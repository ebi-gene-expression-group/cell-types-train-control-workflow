#!/usr/bin/env nextflow 

// initialise required channels 
DATASET_IDS = Channel.create()
NUM_CLUST = Channel.create()
DATA_TYPE = Channel.create()
NORM_METHOD = Channel.create()
BARCODE_COLUMN = Channel.create()
CELL_LABEL_COLUMN = Channel.create()

// parse txt file with dataset accessions into a queue channel 
Channel
    .fromPath(params.data_import.training_dataset_ids)
    .splitCsv(header:false, sep:",")
    .separate(DATASET_IDS, DATA_TYPE, NORM_METHOD, NUM_CLUST, BARCODE_COLUMN, CELL_LABEL_COLUMN)

process fetch_training_datasets {
    publishDir "${baseDir}/data/${dataset_id}", mode: 'copy'
    conda "${baseDir}/envs/load_data.yaml"

    input:
        val(dataset_id) from DATASET_IDS
        val(num_clust) from NUM_CLUST
        val(data_type) from DATA_TYPE
        val(norm_method) from NORM_METHOD

    output:
        set file("data"), val("${dataset_id}") into TRAINING_DATA
        val(num_clust) into N_CLUST

    """
    get_experiment_data.R\
                --accesssion-code ${dataset_id}\
                --config-file ${params.data_import.config_file}\
                --expr-data-type ${data_type}\
                --normalisation-method ${norm_method}\
                --output-dir-name data\
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
        publishDir "${baseDir}/data/${dataset_id}", mode: 'copy'
        conda "${baseDir}/envs/exp_metadata.yaml"

        input:
            set file(data), val(dataset_id) from TRAINING_DATA

        output:
            set file("data"), val("${dataset_id}") into TRAINING_DATA_PROCESSED

        """
        unmelt_condensed.R\
                -i ${data}/condensed-sdrf.tsv\
                -o ${data}/unmelted_sdrf.tsv\
                --has-ontology\
                --retain-types ${params.unmelt_sdrf.retain_types}        
        """
    } 
    COMBINED_TRAINING_DATA = TRAINING_DATA_PROCESSED.merge(BARCODE_COLUMN, CELL_LABEL_COLUMN)
} else {
    COMBINED_TRAINING_DATA = TRAINING_DATA.merge(BARCODE_COLUMN, CELL_LABEL_COLUMN)
}


// fork queue channel contents into channels for corresponding tools
COMBINED_TRAINING_DATA.into{
    GARNETT_TRAINING_DATA
    SCMAP_CELL_TRAINING_DATA
    SCMAP_CLUSTER_TRAINING_DATA
    SCPRED_TRAINING_DATA
}

// add number of clusters to Garnett data 
GARNETT_FULL_DATA = GARNETT_TRAINING_DATA.merge(N_CLUST)

//////////////////////////////////////////
// train each classifier on provided data 
//////////////////////////////////////////

// run garnett training 
if(params.garnett.run == "True"){
    process train_garnett_classifier {
        publishDir "${baseDir}/data/${dataset_id}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }
        
        input:
            set file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col), val(num_clust) from GARNETT_FULL_DATA
            
        output:
            file("garnett_classifier.rds") into GARNETT_CLASSIFIER

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
        """
    }
} else {
    GARNETT_CLASSIFIER = Channel.empty()
}


// run scpred training 
if(params.scpred.run == "True"){
    process run_scpred_workflow {
        publishDir "${baseDir}/data/${dataset_id}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
            set file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col) from SCPRED_TRAINING_DATA

        output:
            file("${params.scpred.trained_model}") into SCPRED_CLASSIFIER

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
        """
    }
} else {
    SCPRED_CLASSIFIER = Channel.empty()
}


// run scmap-cluster training
if(params.scmap_cluster.run == "True"){
    process run_scmap_cluster_workflow {
        publishDir "${baseDir}/data/${dataset_id}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
            set file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col) from SCMAP_CLUSTER_TRAINING_DATA


        output:
            file("scmap_index_cluster.rds")


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
        """
    }
} else {
    SCMAP_CLUSTER_CLASSIFIER = Channel.empty()
}


// run scmap-cell training
if(params.scmap_cell.run == "True"){
    process run_scmap_cell_workflow {
        publishDir "${baseDir}/data/${dataset_id}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
            set file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col) from SCMAP_CELL_TRAINING_DATA


        output:
            file("scmap_index_cell.rds") into SCMAP_CELL_INDEX

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
        """
    }
} else {
    SCMAP_CELL_CLASSIFIER = Channel.empty()
}
