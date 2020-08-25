base = box:Dogon Longitudinal Study/BP and early growth

bptype = sbp

all: models opt_dim growth_patterns_control growth_patterns_nocontrol mediation table1

table1:
	~kshedden/go/bin/rclone copy $(bptype)_table1.csv "$(base)/$(bptype)/mixed"

# use make.jl not this section
mediation:
	pandoc -o $(bptype)/mixed/dim_1/HT_dadbp_med.pdf sbp/mixed/dim_1/HT_dadbp_med.txt
	pandoc -o $(bptype)/mixed/dim_1/HT_dadbp_med.pdf sbp/mixed/dim_1/HT_dadbp_female_med.txt
	pandoc -o $(bptype)/mixed/dim_1/HT_dadbp_med.pdf sbp/mixed/dim_1/HT_dadbp_male_med.txt
	pandoc -o $(bptype)/mixed/dim_1/WT_dadbp_med.pdf sbp/mixed/dim_1/WT_dadbp_med.txt
	pandoc -o $(bptype)/mixed/dim_1/WT_dadbp_med.pdf sbp/mixed/dim_1/WT_dadbp_female_med.txt
	pandoc -o $(bptype)/mixed/dim_1/WT_dadbp_med.pdf sbp/mixed/dim_1/WT_dadbp_male_med.txt
	pandoc -o $(bptype)/mixed/dim_1/BMI_dadbp_med.pdf sbp/mixed/dim_1/BMI_dadbp_med.txt
	pandoc -o $(bptype)/mixed/dim_1/BMI_dadbp_med.pdf sbp/mixed/dim_1/BMI_dadbp_female_med.txt
	pandoc -o $(bptype)/mixed/dim_1/BMI_dadbp_med.pdf sbp/mixed/dim_1/BMI_dadbp_male_med.txt
	rclone copy $(bptype)/mixed/dim_1/HT_dadbp_med.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_1/HT_dadbp_female_med.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_1/HT_dadbp_male_med.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_1/WT_dadbp_med.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_1/WT_dadbp_female_med.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_1/WT_dadbp__male_med.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_1/BMI_dadbp_med.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_1/BMI_dadbp_female_med.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_1/BMI_dadbp_male_med.pdf "$(base)/$(bptype)/mixed/dim_1"

models:
	python bundle.py $(bptype)
	rclone copy $(bptype)/mixed/dim_1/models.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_2/models.pdf "$(base)/$(bptype)/mixed/dim_2"
	python plot_fitted_model.py $(bptype)
	python plot_coeff_traj_pca.py $(bptype) 1
	#python plot_coeff_traj_pca.py $(bptype) 2
	rclone copy $(bptype)/mixed/dim_1/model_plots.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_1/coeff_traj_pca.pdf "$(base)/$(bptype)/mixed/dim_1"
	rclone copy $(bptype)/mixed/dim_2/coeff_traj_pca.pdf "$(base)/$(bptype)/mixed/dim_2"
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

growth_patterns_nocontrol:
	python growth_patterns_nocontrol.py $(bptype) 1 3 nodadbp
	python growth_patterns_nocontrol.py $(bptype) 1 3 dadbp
	python growth_patterns_nocontrol_plot.py $(bptype) 1 nodadbp
	python growth_patterns_nocontrol_plot.py $(bptype) 1 dadbp
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_nocontrolcbs_nodadbp_1_3.pdf $(bptype)/mixed/dim_1/growth_patterns_nocontrolcbs_nodadbp_1_3.txt
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_nocontrolcbs_dadbp_1_3.pdf $(bptype)/mixed/dim_1/growth_patterns_nocontrolcbs_dadbp_1_3.txt
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_nocontrolcbs_nodadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_nocontrolcbs_dadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_plot_nocontrolcbs_nodadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_plot_nocontrolcbs_dadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"

growth_patterns_control:
	python growth_patterns.py $(bptype) 1 1 none nodadbp
	python growth_patterns.py $(bptype) 1 3 none nodadbp
	python growth_patterns.py $(bptype) 1 3 none dadbp
	python growth_patterns.py $(bptype) 1 4 none nodadbp
	python growth_patterns.py $(bptype) 1 3 zf nodadbp
	python growth_patterns.py $(bptype) 1 3 zm nodadbp
	python growth_patterns_plot.py $(bptype) 1 nodadbp
	python growth_patterns_plot.py $(bptype) 1 dadbp
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_1.pdf $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_1.txt
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_3_zf.pdf $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_3_zf.txt
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_3_zm.pdf $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_3_zm.txt
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_3.pdf $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_3.txt
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_controlcbs_dadbp_1_3.pdf $(bptype)/mixed/dim_1/growth_patterns_controlcbs_dadbp_1_3.txt
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_4.pdf $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_4.txt
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_controlcbs_dadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_4.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_3_zm.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_1.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_3_zf.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_plot_controlcbs_nodadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
	~kshedden/go/bin/rclone copy $(bptype)/mixed/dim_1/growth_patterns_plot_controlcbs_dadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
