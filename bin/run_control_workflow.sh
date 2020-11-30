#!/usr/bin/env bash 

export WORKFLOW_ROOT=`pwd`
#conda activate nextflow
if [ $2 = 'resume' ]; then
    resume='-resume'
else
    resume=''
fi
 
nextflow run main.nf -profile $1 $resume
