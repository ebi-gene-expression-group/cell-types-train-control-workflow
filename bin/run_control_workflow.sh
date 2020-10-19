#!/usr/bin/env bash 

export WORKFLOW_ROOT=`pwd`
#conda activate nextflow
nextflow run main.nf -profile $1 -resume
