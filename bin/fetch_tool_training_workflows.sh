#!/usr/bin/env bash
TRAIN_WORKFLOW_ROOT="$PWD"
TRAIN_WORKFLOWS="$PWD/cell-types-train-workflows"

# Clone or update the cell-types-train-workflows repo containing submodules for individual pipelines
if [ ! -d 'cell-types-train-workflows' ]; then
    git clone --recursive https://github.com/ebi-gene-expression-group/cell-types-train-workflows.git $TRAIN_WORKFLOWS
fi

pushd $TRAIN_WORKFLOWS > /dev/null
git checkout origin/develop > /dev/null
git pull origin develop > /dev/null
git submodule update --recursive --remote > /dev/null
popd > /dev/null

