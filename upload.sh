
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

pandoc -o $1/$2/opt_dim_0.pdf $1/$2/opt_dim_0.txt
pandoc -o $1/$2/opt_dim_1.pdf $1/$2/opt_dim_1.txt

python bundle.py $1 $2

base="remote:Dogon Longitudinal Study/SBP and early growth"

python plot_fitted_model.py $1 $2
rclone copy model_plots.pdf "${base}/$1/$2/dim_1"

rclone copy $1/$2/dim_1/coeff_traj_pca.pdf "${base}/$1/$2/dim_1"
rclone copy $1/$2/dim_2/coeff_traj_pca.pdf "${base}/$1/$2/dim_2"
##rclone copy $1/$2/dim_3/coeff_traj_pca.pdf "${base}/$1/$2/dim_3"

rclone copy $1/$2/dim_1/models_0.pdf "${base}/$1/$2/dim_1"
rclone copy $1/$2/dim_1/models_1.pdf "${base}/$1/$2/dim_1"

rclone copy $1/$2/dim_2/models_0.pdf "${base}/$1/$2/dim_2"
rclone copy $1/$2/dim_2/models_1.pdf "${base}/$1/$2/dim_2"

rclone copy $1/$2/dim_3/models_0.pdf "${base}/$1/$2/dim_3"
rclone copy $1/$2/dim_3/models_1.pdf "${base}/$1/$2/dim_3"

rclone copy $1/$2/opt_dim_0.pdf "${base}/$1/$2"
rclone copy $1/$2/opt_dim_1.pdf "${base}/$1/$2"

python merge_stats.py $1 $2
rclone copy sample_sizes.csv "${base}/$1/$2"
