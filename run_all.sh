#!/bin/bash

if [ -z "$1" ]
then
    echo "Required arguments: [sbp|dbp] [gee|mixed]"
    exit 1
fi

if [ -z "$2" ]
then
    echo "Required arguments: [sbp|dbp] [gee|mixed]"
    exit 1
fi

bash sbp_pca.sh $1 $2 1 growth >&t1
bash sbp_pca.sh $1 $2 1 nogrowth >&t2
python plot_coeff_traj_pca.py $1 $2 1 >&t3

bash sbp_pca.sh $1 $2 2 growth >&t6
bash sbp_pca.sh $1 $2 2 nogrowth >&t7
python plot_coeff_traj_pca.py $1 $2 2 >&t8

