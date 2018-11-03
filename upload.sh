
python opt_dim.py $1
python bundle.py $1

pandoc -o $1/opt_dim_0.pdf $1/opt_dim_0.txt
pandoc -o $1/opt_dim_1.pdf $1/opt_dim_1.txt

python plot_fitted_model.py mixed

rclone copy $1/dim_1/coeff_traj_pca.pdf  "remote:Dogon Longitudinal Study/SBP and early growth/$1/dim_1"
rclone copy $1/dim_2/coeff_traj_pca.pdf  "remote:Dogon Longitudinal Study/SBP and early growth/$1/dim_2"
rclone copy $1/dim_3/coeff_traj_pca.pdf  "remote:Dogon Longitudinal Study/SBP and early growth/$1/dim_3"

rclone copy $1/dim_1/models_0.pdf  "remote:Dogon Longitudinal Study/SBP and early growth/$1/dim_1"
rclone copy $1/dim_1/models_1.pdf  "remote:Dogon Longitudinal Study/SBP and early growth/$1/dim_1"

rclone copy $1/dim_2/models_0.pdf  "remote:Dogon Longitudinal Study/SBP and early growth/$1/dim_2"
rclone copy $1/dim_2/models_1.pdf  "remote:Dogon Longitudinal Study/SBP and early growth/$1/dim_2"

rclone copy $1/dim_3/models_0.pdf  "remote:Dogon Longitudinal Study/SBP and early growth/$1/dim_3"
rclone copy $1/dim_3/models_1.pdf  "remote:Dogon Longitudinal Study/SBP and early growth/$1/dim_3"

rclone copy $1/opt_dim_0.pdf "remote:Dogon Longitudinal Study/SBP and early growth/$1"
rclone copy $1/opt_dim_1.pdf "remote:Dogon Longitudinal Study/SBP and early growth/$1"

python merge_stats.py
rclone copy sample_sizes.csv "remote:Dogon Longitudinal Study/SBP and early growth/$1"
