#!/usr/bin/env nextflow 

// parse txt file with dataset accessions into a queue channel 
DATASET_IDS = Channel.create()
NUM_CLUST = Channel.create()
BARCODE_COLUMN = Channel.create()
CELL_LABEL_COLUMN = Channel.create()

Channel
    .fromPath(params.training_dataset_ids)
    .splitCsv(header:false, sep:",")
    .separate(DATASET_IDS, NUM_CLUST, BARCODE_COLUMN, CELL_LABEL_COLUMN)

process fetch_training_datasets {
    conda "${baseDir}/envs/load_data.yaml"

    input:
        val(dataset_id) from DATASET_IDS
        val(num_clust) from NUM_CLUST

    output:
        set file("${dataset_id}"), val("${num_clust}") into TRAINING_DATASET

    """
    get_experiment_data.R\
                --accesssion-code ${dataset_id}\
                --config-file ${params.data_import.config_file}\
                --expr-data-type ${params.data_import.expr_data_type}\
                --normalisation-method ${params.data_import.norm_method}\
                --output-dir-name ${dataset_id}\
                --get-sdrf ${params.data_import.get_sdrf}\
                --get-condensed-sdrf ${params.data_import.get_cond_sdrf}\
                --get-idf ${params.data_import.get_idf}\
                --get-marker-genes ${params.data_import.get_marker_genes}\
                --number-of-clusters ${num_clust}
    """
}

COMBINED_TRAINING_DATA = TRAINING_DATASET.merge(BARCODE_COLUMN, CELL_LABEL_COLUMN)
//COMBINED_TRAINING_DATA.view()



// duplicate queue channel contents into corresponding tool channels 



//////////////////////////////////////////
// train each classifier on provided data 
//////////////////////////////////////////


// run garnett training 
if(params.garnett.run == "True"){
    process train_garnett_classifier {
        publishDir "${}", mode: 'copy' //TODO: specify piblishDir 
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }
        
        input:
            file(training_data) from GARNETT_TRAINING_DATA
            val(training_id) from GARNETT_TRAINING_IDS
        
        output:
            file("garnett_classifier.rds") into GARNETT_CLASSIFIER

        """
        RESULTS_DIR=\$PWD

        nextflow run $TRAIN_WORKFLOWS/garnett-train-workflow/main.nf\
                            -profile cluster\
                            --results_dir \$RESULTS_DIR\
                            --ref_10x_dir ${reference_10X_dir}\
                            --marker_genes ${ref_markers}\
                            --marker_gene_id_type ${params.garnett.marker_gene_id_type}\
                            --classifier_gene_type ${params.garnett.classifier_gene_type}\
                            --n_outgroups ${params.garnett.n_outgroups}\
        """
    }
}else{
    GARNETT_CLASSIFIER = Channel.empty()
}


// run scpred training 
// run scpred 
if(params.scpred.run == "True"){
    process run_scpred_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:



        output:

        """
        """
    }
}else{
    SCPRED_CLASSIFIER = Channel.empty()
}











