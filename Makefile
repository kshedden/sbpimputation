
base = remote:Dogon Longitudinal Study/BP and early growth

bptype = sbp

all: models opt_dim growth_patterns

table1:
	~kshedden/go/bin/rclone copy $(bptype)_table1.csv "$(base)/$(bptype)/mixed"

models:
	python bundle.py $(bptype) mixed
	python plot_fitted_model.py $(bptype) mixed
	rclone copy $(bptype)/mixed/dim_1/model_plots.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_1/coeff_traj_pca.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_2/coeff_traj_pca.pdf "$(base)/$(bptype)/mixed/dim_2"
	rclone copy $(bptype)/mixed/dim_1/models.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_2/models.pdf "$(base)/$(bptype)/mixed/dim_2"
	rclone copy $(bptype)/mixed/stats/HT_n.csv "$(base)/$(bptype)/mixed/stats"
	rclone copy $(bptype)/mixed/stats/WT_n.csv "$(base)/$(bptype)/mixed/stats"
	rclone copy $(bptype)/mixed/stats/BMI_n.csv "$(base)/$(bptype)/mixed/stats"
	python excluded_ids.py
	rclone copy excluded_ids.csv "$(base)/$(bptype)/mixed/stats"
	rclone copy ht_stats.csv "$(base)/$(bptype)/mixed/stats"
	rclone copy wt_stats.csv "$(base)/$(bptype)/mixed/stats"
	rclone copy bmi_stats.csv "$(base)/$(bptype)/mixed/stats"

opt_dim:
	python opt_dim.py $(bptype) mixed
	pandoc -V geometry:landscape -o $(bptype)/mixed/opt_dim.pdf $(bptype)/mixed/opt_dim.txt
	rclone copy $(bptype)/mixed/opt_dim.pdf "$(base)/$(bptype)/mixed"

growth_patterns:
	python growth_patterns.py $(bptype) 1 1 none
	python growth_patterns.py $(bptype) 1 3 none
	python growth_patterns.py $(bptype) 1 4 none
	python growth_patterns.py $(bptype) 1 3 zf
	python growth_patterns.py $(bptype) 1 3 zm
	python growth_patterns_plot.py $(bptype) 1
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_1_1.pdf $(bptype)/mixed/dim_1/growth_patterns_1_1.txt
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_1_3_zf.pdf $(bptype)/mixed/dim_1/growth_patterns_1_3_zf.txt
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_1_3_zm.pdf $(bptype)/mixed/dim_1/growth_patterns_1_3_zm.txt
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_1_3.pdf $(bptype)/mixed/dim_1/growth_patterns_1_3.txt
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_1_4.pdf $(bptype)/mixed/dim_1/growth_patterns_1_4.txt
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_1_1.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_1_3_zf.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_1_3_zm.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_1_4.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_plot_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
