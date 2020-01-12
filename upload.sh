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

python opt_dim.py $1 $2

pandoc -o $1/$2/opt_dim.pdf $1/$2/opt_dim.txt

python bundle.py $1 $2

base="remote:Dogon Longitudinal Study/BP and early growth"

python plot_fitted_model.py $1 $2
rclone copy $1/$2/dim_1/model_plots.pdf "${base}/$1/$2/dim_1"

rclone copy $1/$2/dim_1/coeff_traj_pca.pdf "${base}/$1/$2/dim_1"
rclone copy $1/$2/dim_2/coeff_traj_pca.pdf "${base}/$1/$2/dim_2"

rclone copy $1/$2/dim_1/models.pdf "${base}/$1/$2/dim_1"
rclone copy $1/$2/dim_2/models.pdf "${base}/$1/$2/dim_2"

rclone copy $1/$2/opt_dim.pdf "${base}/$1/$2"

python merge_stats.py $1 $2
rclone copy sample_sizes.csv "${base}/$1/$2"

python growth_patterns.py $1 1
pandoc -o $1/$2/dim_1/growth_patterns.pdf $1/$2/dim_1/growth_patterns.txt
rclone copy $1/$2/dim_1/growth_patterns.pdf "${base}/$1/$2/dim_1"

python growth_patterns_plot.py $1 1
rclone copy $1/$2/dim_1/growth_patterns_plot.pdf "${base}/$1/$2/dim_1"
