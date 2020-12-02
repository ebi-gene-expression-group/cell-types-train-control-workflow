#!/usr/bin/env nextflow 

// extract matrix types required by the tools, conditioned on them being 'on' 
tool_switch = ["True":0, "False":1]
garnett_matrix_type = [params.garnett.matrix_type, null][tool_switch[params.garnett.run]]
scmap_cluster_matrix_type = [params.scmap_cluster.matrix_type, null][tool_switch[params.scmap_cluster.run]]
scmap_cell_matrix_type = [params.scmap_cell.matrix_type, null][tool_switch[params.scmap_cell.run]]
scpred_matrix_type = [params.scpred.matrix_type, null][tool_switch[params.scpred.run]]

UNIQUE_MATRIX_TYPES = Channel
                    .of(garnett_matrix_type,
                        scmap_cluster_matrix_type,
                        scmap_cell_matrix_type,
                        scpred_matrix_type)
                    .filter{ it != null }
                    .unique()

// parse txt file with dataset accessions into a queue channel; build combinations with matrix types 
IMPORT_PARAMS = Channel
                .fromPath(params.data_import.training_datasets)
                .splitCsv(header:false, sep:",")
                .combine(UNIQUE_MATRIX_TYPES)

process fetch_training_datasets {
    conda "${baseDir}/envs/load_data.yaml"

    input:
        tuple val(dataset_id), val(seq_method), val(num_clust), val(barcode_col), val(cell_label_col), val(matrix_type) from IMPORT_PARAMS

    output:
        tuple file(dataset_id), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) into TRAINING_DATA
        val(num_clust) into N_CLUST

    """
    if [ ${seq_method} ==  "droplet" ]; then 
        MATRIX_TYPE_UPD="CPM"
    else
        MATRIX_TYPE_UPD=${matrix_type}
    fi

    get_experiment_data.R\
                --accesssion-code ${dataset_id}\
                --config-file ${params.data_import.scxa_import_config_file}\
                --matrix-type \$MATRIX_TYPE_UPD\
                --output-dir-name ${dataset_id}\
                --get-sdrf ${params.data_import.get_sdrf}\
                --get-condensed-sdrf ${params.data_import.get_cond_sdrf}\
                --get-idf ${params.data_import.get_idf}\
                --get-marker-genes ${params.data_import.get_marker_genes}\
                --number-of-clusters ${num_clust}
    """
}

// filter out unwanted characters from SDRF and marker gene files
process filter_labels{
    conda "${baseDir}/envs/label_analysis.yaml"

    input: 
        tuple file(data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) from TRAINING_DATA
        val(num_clust) from N_CLUST

    output:
        tuple file("${data}_upd"), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) into TRAINING_DATA_FILT
        val(num_clust) into N_CLUST_FILT

    """
    # SDRF file
    check_labels.R\
          --input-file ${data}/condensed-sdrf.tsv\
          --label-field '${cell_label_col}'\
          --condensed\
          --output-path ${data}/condensed_sdrf_filt.tsv &&
    check_labels.R\
          --input-file ${data}/marker_genes_${num_clust}.tsv\
          --label-field '${params.garnett.groups_col}'\
          --output-path ${data}/marker_genes_${num_clust}_filt.tsv
    
    mv ${data} "${data}_upd"
    """

}


// to avoid problems with sdrf-barcode matching, unmelt condensed SDRF and use it in downstream processes
if(params.unmelt_sdrf.run == "True"){
    process unmelt_condensed_sdrf {
        conda "${baseDir}/envs/exp_metadata.yaml"

        memory { 10.GB * task.attempt }
        maxRetries 5
        errorStrategy { task.attempt<=5 ? 'retry' : 'ignore' }

        input:
            tuple file(data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) from TRAINING_DATA_FILT

        output:
            tuple file(data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) into TRAINING_DATA_UNMELT

        """
        unmelt_condensed.R\
                -i ${data}/condensed_sdrf_filt.tsv\
                -o ${data}/unmelted_sdrf.tsv\
                --has-ontology ${params.unmelt_sdrf.has_ontology}\
                --retain-types ${params.unmelt_sdrf.retain_types}        
        """
    } 
    TRAINING_DATA_UNMELT.set{ TRAINING_DATA_PROCESSED }
} else {
    TRAINING_DATA.set{ TRAINING_DATA_PROCESSED }
}

// fork queue channel contents into channels for corresponding tools
TRAINING_DATA_PROCESSED.into{
    GARNETT_TRAINING_DATA
    SCMAP_CELL_TRAINING_DATA
    SCMAP_CLUSTER_TRAINING_DATA
    SCPRED_TRAINING_DATA
}

// add number of clusters to Garnett data 
GARNETT_FULL_DATA = N_CLUST_FILT.merge(GARNETT_TRAINING_DATA)

//////////////////////////////////////////
// train each classifier on provided data 
//////////////////////////////////////////

// keep only relevant version of the dataset 
GARNETT_FILTERED_DATA = GARNETT_FULL_DATA.filter{ it[5] == params.garnett.matrix_type }

// run garnett training 
if(params.garnett.run == "True"){
    process run_garnett_workflow {
        publishDir "${baseDir}/data/${dataset_id}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.attempt<=3  ? 'retry' : 'ignore' }   
        maxRetries 3
        memory { 16.GB * task.attempt }
     
        maxForks 5

        input:
            tuple val(num_clust), file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) from GARNETT_FILTERED_DATA

        output:
            file("${dataset_id}_garnett.rds") into GARNETT_CLASSIFIER

        """
        RESULTS_DIR=\$PWD

        nextflow run $TRAIN_WORKFLOWS/garnett-train-workflow/main.nf\
			                -profile ${params.profile}\
                            --results_dir \$RESULTS_DIR\
                            --training_10x_dir ${training_data}/10x_data\
                            --training_dataset_id ${dataset_id}\
                            --marker_genes ${training_data}/marker_genes_${num_clust}_filt.tsv\
                            --pval_col ${params.garnett.pval_col}\
                            --groups_col ${params.garnett.groups_col}\
                            --gene_names ${params.garnett.gene_names}\
                            --database ${params.garnett.database}\
                            --marker_gene_id_type ${params.garnett.marker_gene_id_type}\
                            --classifier_gene_type ${params.garnett.classifier_gene_type}\
                            --n_outgroups ${params.garnett.n_outgroups}

        mv garnett_classifier.rds ${dataset_id}_garnett.rds
        """
    }
} else {
    GARNETT_CLASSIFIER = Channel.empty()
}


// keep only relevant version of the dataset 
SCPRED_FILTERED_DATA = SCPRED_TRAINING_DATA.filter{ it[4] == params.scpred.matrix_type }
// run scpred training
if(params.scpred.run == "True"){
    process run_scpred_workflow {
        publishDir "${baseDir}/data/${dataset_id}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.attempt<=3  ? 'retry' : 'ignore' }
        maxRetries 3
        memory { 16.GB * task.attempt }

        maxForks 5
        
        input:
            tuple file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) from SCPRED_FILTERED_DATA

        output:
            file("${dataset_id}_scpred.rds") into SCPRED_CLASSIFIER

        """
        RESULTS_DIR=\$PWD

        nextflow run $TRAIN_WORKFLOWS/scpred-train-workflow/main.nf\
                            -profile ${params.profile}\
			                --results_dir \$RESULTS_DIR\
                            --exclusions ${params.exclusions}\
                            --method ${params.scpred.method}\
                            --training_10x_dir ${training_data}/10x_data\
                            --metadata_file ${training_data}/unmelted_sdrf.tsv\
                            --training_dataset_id ${dataset_id}\
                            --train_probs_plot_path ${params.scpred.train_probs_plot_path}\
                            --normalised_counts_slot ${params.scpred.normalised_counts_slot}\
                            --cell_id_col_name "${barcode_col}"\
                            --cell_types_col_name "${cell_label_col}"\
                            --col_names ${params.scpred.col_names}\
                            --trained_model ${params.scpred.trained_model}\
                            --log_transform ${params.scpred.log_transform}\
                            --model ${params.scpred.model}

        mv scpred_classifier.rds ${dataset_id}_scpred.rds
        """
    }
} else {
    SCPRED_CLASSIFIER = Channel.empty()
}


// keep only relevant version of the dataset 
SCMAP_CLUSTER_FILTERED_DATA = SCMAP_CLUSTER_TRAINING_DATA.filter{ it[4] == params.scmap_cluster.matrix_type }
// run scmap-cluster training
if(params.scmap_cluster.run == "True"){
    process run_scmap_cluster_workflow {
        publishDir "${baseDir}/data/${dataset_id}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.attempt<=3  ? 'retry' : 'ignore' }
        maxRetries 3
        memory { 16.GB * task.attempt }

        maxForks 5
        
        input:
            tuple file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) from SCMAP_CLUSTER_FILTERED_DATA

        output:
            file("${dataset_id}_scmap-cluster.rds")

        """
        RESULTS_DIR=\$PWD

        nextflow run $TRAIN_WORKFLOWS/scmap-train-workflow/main.nf\
			                -profile ${params.profile}\
                            --results_dir \$RESULTS_DIR\
                            --training_10x_dir ${training_data}/10x_data\
                            --training_metadata ${training_data}/unmelted_sdrf.tsv\
                            --projection_method cluster\
                            --training_dataset_id ${dataset_id}\
                            --col_names ${params.scmap_cluster.col_names}\
                            --exclusions ${params.exclusions}\
                            --cell_id_col ${barcode_col}\
                            --cluster_col ${cell_label_col}\
                            --threshold ${params.scmap_cluster.threshold}

        mv scmap_index_cluster.rds ${dataset_id}_scmap-cluster.rds
        """
    }
} else {
    SCMAP_CLUSTER_CLASSIFIER = Channel.empty()
}

// keep only relevant version of the dataset 
SCMAP_CELL_FILTERED_DATA = SCMAP_CELL_TRAINING_DATA.filter{ it[4] == params.scmap_cell.matrix_type }
// run scmap-cell training
if(params.scmap_cell.run == "True"){
    process run_scmap_cell_workflow {
        publishDir "${baseDir}/data/${dataset_id}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.attempt<=3  ? 'retry' : 'ignore' }
        maxRetries 3
        memory { 16.GB * task.attempt }

        maxForks 5

        input:
            tuple file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) from SCMAP_CELL_FILTERED_DATA


        output:
            file("${dataset_id}_scmap-cell.rds") into SCMAP_CELL_INDEX

        """
        RESULTS_DIR=\$PWD

        nextflow run $TRAIN_WORKFLOWS/scmap-train-workflow/main.nf\
                            -profile ${params.profile}\
			                --results_dir \$RESULTS_DIR\
                            --training_10x_dir ${training_data}/10x_data\
                            --training_metadata ${training_data}/unmelted_sdrf.tsv\
                            --projection_method cell\
                            --training_dataset_id ${dataset_id}\
                            --col_names ${params.scmap_cell.col_names}\
                            --cell_id_col ${barcode_col}\
                            --cluster_col ${cell_label_col}\
                            --exclusions ${params.exclusions}\
                            --threshold ${params.scmap_cell.threshold}

        mv scmap_index_cell.rds ${dataset_id}_scmap-cell.rds
        """
    }
} else {
    SCMAP_CELL_CLASSIFIER = Channel.empty()
}
