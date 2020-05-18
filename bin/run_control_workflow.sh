#!/usr/bin/env bash 

export WORKFLOW_ROOT=`pwd`

nextflow run main.nf -profile cluster
