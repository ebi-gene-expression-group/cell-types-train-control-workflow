# cell-types-train-control-workflow
Schematic representation of the process is shown below: 
![](classifier_training.png)

This workflow scans a comma-separated config file for specified SCXA dataest accession numbers, imports them and trains a range of classifiers for each dataset. The following columns are expected in the config file: 
* `dataset id`
* `technology type` _("droplet" or "smart-seq")_
* `matrix type` _(raw, filtered, CPM- or TPM-normalised)_
* `number of clusters in marker gene file` 
* `barcode column` (in SDRF file) 
* `cell type column` 

### Running the workflow 
Prior to running the workflow, you will need to fetch and update the submodules for individual tool workflows. Run the following from the workflow directory: `./bin/fetch_tool_training_workflows.sh`. 

You will need Bioconda installed to run the workflow. It is recommended to use a clean environment to run the process. Issue the following commands: 
```
conda create -n nextflow && conda activate nextflow 
conda install nextflow 
./bin/run_control_workflow.sh 
```
In the `run_control_workflow.sh` script, the `<profile>` parameter might be either `local` or `cluster` depending on where you run the process. 
