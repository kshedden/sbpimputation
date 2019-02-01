

base = remote:Dogon Longitudinal Study/BP and early growth

sbp_models:
	python plot_fitted_model.py sbp mixed
	rclone copy sbp/mixed/dim_1/model_plots.pdf "$(base)/sbp/mixed/dim_1"
	rclone copy sbp/mixed/dim_1/coeff_traj_pca.pdf "$(base)/sbp/mixed/dim_1"
	rclone copy sbp/mixed/dim_2/coeff_traj_pca.pdf "$(base)/sbp/mixed/dim_2"
	rclone copy sbp/mixed/dim_1/models.pdf "$(base)/sbp/mixed/dim_1"
	rclone copy sbp/mixed/dim_2/models.pdf "$(base)/sbp/mixed/dim_2"

sbp_opt_dim:
	python opt_dim.py sbp mixed
	pandoc -V geometry:landscape -o sbp/mixed/opt_dim.pdf sbp/mixed/opt_dim.txt
	rclone copy sbp/mixed/opt_dim.pdf "$(base)/sbp/mixed"

dbp_opt_dim:
	python opt_dim.py dbp mixed
	pandoc -V geometry:landscape -o dbp/mixed/opt_dim.pdf dbp/mixed/opt_dim.txt
	rclone copy dbp/mixed/opt_dim.pdf "$(base)/dbp/mixed"

sbp_growth_patterns:
	python growth_patterns.py sbp 1
	python growth_patterns_plot.py sbp 1
	pandoc -o sbp/mixed/dim_1/growth_patterns.pdf sbp/mixed/dim_1/growth_patterns.txt
	rclone copy sbp/mixed/dim_1/growth_patterns.pdf "$(base)/sbp/mixed/dim_1"
	rclone copy sbp/mixed/dim_1/growth_patterns_plot.pdf "$(base)/sbp/mixed/dim_1"

dbp_growth_patterns:
	python growth_patterns.py dbp 1
	python growth_patterns_plot.py dbp 1
	pandoc -o dbp/mixed/dim_1/growth_patterns.pdf dbp/mixed/dim_1/growth_patterns.txt
	rclone copy dbp/mixed/dim_1/growth_patterns.pdf "${base}/dbp/mixed/dim_1"
	rclone copy dbp/mixed/dim_1/growth_patterns_plot.pdf "${base}/dbp/mixed/dim_1"
