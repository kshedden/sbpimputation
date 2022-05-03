#!/bin/bash

if [ -z "$1" ]
then
    echo "Required arguments: [sbp|dbp]"
    exit 1
fi

python run_sbp_pca_sex.py $1 1 >& t1
