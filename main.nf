#!/usr/bin/env nextflow 

// initialise required channels 
DATASET_IDS = Channel.create()
NUM_CLUST = Channel.create()
BARCODE_COLUMN = Channel.create()
CELL_LABEL_COLUMN = Channel.create()

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
        tuple val(dataset_id), val(seq_method), val(num_clust), val(barcode_col), val(cell_type_col), val(matrix_type) from IMPORT_PARAMS

    output:
        tuple file("data"), val(dataset_id), val(barcode_col), val(cell_type_col), val(matrix_type) into TRAINING_DATA
        val(num_clust) into N_CLUST

    """
    if [ ${seq_method} ==  "droplet" ]; then 
        MATRIX_TYPE_UPD="CPM"
    else
        MATRIX_TYPE_UPD=${matrix_type}
    fi

    get_experiment_data.R\
                --accesssion-code ${dataset_id}\
                --config-file ${params.data_import.config_file}\
                --matrix-type \$MATRIX_TYPE_UPD\
                --output-dir-name ${dataset_id}\
                --get-sdrf ${params.data_import.get_sdrf}\
                --get-condensed-sdrf ${params.data_import.get_cond_sdrf}\
                --get-idf ${params.data_import.get_idf}\
                --get-marker-genes ${params.data_import.get_marker_genes}\
                --number-of-clusters ${num_clust}
    """
}

// to avoid problems with sdrf-barcode matching, unmelt condensed SDRF and use it in downstream processes
if(params.unmelt_sdrf.run == "True"){
    process unmelt_condensed_sdrf {
        conda "${baseDir}/envs/exp_metadata.yaml"

        input:
            tuple file(data), val(dataset_id), val(barcode_col), val(cell_type_col), val(matrix_type) from TRAINING_DATA

        output:
            tuple file(data), val(dataset_id), val(barcode_col), val(cell_type_col), val(matrix_type) into TRAINING_DATA_UNMELT

        """
        unmelt_condensed.R\
                -i ${data}/condensed-sdrf.tsv\
                -o ${data}/unmelted_sdrf.tsv\
                --has-ontology\
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
GARNETT_FULL_DATA = N_CLUST.merge(GARNETT_TRAINING_DATA)

//////////////////////////////////////////
// train each classifier on provided data 
//////////////////////////////////////////

// keep only relevant version of the dataset 
GARNETT_FILTERED_DATA = GARNETT_FULL_DATA.filter{ it[5] == params.garnett.matrix_type }

// run garnett training 
if(params.garnett.run == "True"){
    process train_garnett_classifier {
        publishDir "${baseDir}/data/${dataset_id}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }
        
        input:
            tuple val(num_clust), file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) from GARNETT_FILTERED_DATA
            
        output:
            file("garnett_classifier_${dataset_id}.rds") into GARNETT_CLASSIFIER

        """
        RESULTS_DIR=\$PWD

        nextflow run $TRAIN_WORKFLOWS/garnett-train-workflow/main.nf\
                            -profile ${params.profile}\
                            --results_dir \$RESULTS_DIR\
                            --training_10x_dir ${training_data}/10x_data\
                            --training_dataset_id ${dataset_id}\
                            --marker_genes ${training_data}/marker_genes_${num_clust}.tsv\
                            --pval_col ${params.garnett.pval_col}\
                            --groups_col ${params.garnett.groups_col}\
                            --gene_names ${params.garnett.gene_names}\
                            --database ${params.garnett.database}\
                            --marker_gene_id_type ${params.garnett.marker_gene_id_type}\
                            --classifier_gene_type ${params.garnett.classifier_gene_type}\
                            --n_outgroups ${params.garnett.n_outgroups}

        mv garnett_classifier.rds garnett_classifier_${dataset_id}.rds
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

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
            tuple file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) from SCPRED_FILTERED_DATA

        output:
            file("scpred_classifier_${dataset_id}.rds") into SCPRED_CLASSIFIER

        """
        RESULTS_DIR=\$PWD

        nextflow run $TRAIN_WORKFLOWS/scpred-train-workflow/main.nf\
                            -profile ${params.profile}\
                            --results_dir \$RESULTS_DIR\
                            --method ${params.scpred.method}\
                            --training_10x_dir ${training_data}/10x_data\
                            --metadata_file ${training_data}/unmelted_sdrf.tsv\
                            --training_dataset_id ${dataset_id}\
                            --eigenvalue_plot_path ${params.scpred.eigenvalue_plot_path}\
                            --train_probs_plot_path ${params.scpred.train_probs_plot_path}\
                            --normalised_counts_slot ${params.scpred.normalised_counts_slot}\
                            --cell_id_col_name ${barcode_col}\
                            --cell_types_col_name ${cell_label_col}\
                            --col_names ${params.scpred.col_names}\
                            --trained_model ${params.scpred.trained_model}\
                            --log_transform ${params.scpred.log_transform}\
                            --model ${params.scpred.model}

        mv scpred_classifier.rds scpred_classifier_${dataset_id}.rds
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

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
            tuple file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) from SCMAP_CLUSTER_FILTERED_DATA

        output:
            file("scmap_index_cluster_${dataset_id}.rds")

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
                            --cell_id_col ${barcode_col}\
                            --cluster_col ${cell_label_col}\
                            --threshold ${params.scmap_cluster.threshold}

        mv scmap_index_cluster.rds scmap_index_cluster_${dataset_id}.rds
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

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
            tuple file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) from SCMAP_CELL_FILTERED_DATA


        output:
            file("scmap_index_cell_${dataset_id}.rds") into SCMAP_CELL_INDEX

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
                            --threshold ${params.scmap_cell.threshold}

        mv scmap_index_cell.rds scmap_index_cell_${dataset_id}.rds
        """
    }
} else {
    SCMAP_CELL_CLASSIFIER = Channel.empty()
}
