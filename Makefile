base=Dogon Longitudinal Study/BP_and_early_growth2

DBU=dropbox_uploader.sh

bptype = sbp

all: models opt_dim growth_patterns_control growth_patterns_nocontrol mediation table1

growth_education:
	$(DBU) upload growth_educ.txt "$(base)"
	$(DBU) upload growth_loadings.pdf "$(base)"
	$(DBU) upload growth_educ_cox.txt "$(base)"
	$(DBU) upload growth_loadings_cox.pdf "$(base)"

tracking:
	$(DBU) upload tracking.pdf "$(base)"

summary:
	$(DBU) upload summary.txt "$(base)/$(bptype)/mixed"
	$(DBU) upload ref_ids_female.csv "$(base)/$(bptype)/mixed"
	$(DBU) upload ref_ids_male.csv "$(base)/$(bptype)/mixed"
	$(DBU) upload haz_histogram.pdf "$(base)"

table1:
	$(DBU) upload $(bptype)_table1.csv "$(base)/$(bptype)/mixed"

cohort_plots:
	$(DBU) upload cohort_plots.pdf "$(base)"

pca_basis:
	$(DBU) upload $(bptype)/mixed/dim_1/HT_basis.pdf "$(base)/$(bptype)/mixed/dim_1"
	$(DBU) upload $(bptype)/mixed/dim_1/WT_basis.pdf "$(base)/$(bptype)/mixed/dim_1"
	$(DBU) upload $(bptype)/mixed/dim_1/BMI_basis.pdf "$(base)/$(bptype)/mixed/dim_1"

models:
	python bundle.py $(bptype)
	$(DBU) upload $(bptype)/mixed/dim_1/models.txt "$(base)/$(bptype)/mixed/dim_1"
	#$(DBU) upload $(bptype)/mixed/dim_2/models.txt "$(base)/$(bptype)/mixed/dim_2"
	#python plot_fitted_model.py $(bptype)
	#python plot_coeff_traj_pca.py $(bptype) 1
	#python plot_coeff_traj_pca.py $(bptype) 2
	#$(DBU) upload $(bptype)/mixed/dim_1/model_plots.pdf "$(base)/$(bptype)/mixed/dim_1"
	#$(DBU) upload $(bptype)/mixed/dim_1/coeff_traj_pca.pdf "$(base)/$(bptype)/mixed/dim_1"
	#$(DBU) upload $(bptype)/mixed/dim_2/coeff_traj_pca.pdf "$(base)/$(bptype)/mixed/dim_2"
	#$(DBU) upload $(bptype)/mixed/stats/HT_n.csv "$(base)/$(bptype)/mixed/stats"
	#$(DBU) upload $(bptype)/mixed/stats/WT_n.csv "$(base)/$(bptype)/mixed/stats"
	#$(DBU) upload $(bptype)/mixed/stats/BMI_n.csv "$(base)/$(bptype)/mixed/stats"
	#python excluded_ids.py
	#$(DBU) upload excluded_ids.csv "$(base)/$(bptype)/mixed/stats"
	#$(DBU) upload ht_stats.csv "$(base)/$(bptype)/mixed/stats"
	#$(DBU) upload wt_stats.csv "$(base)/$(bptype)/mixed/stats"
	#$(DBU) upload bmi_stats.csv "$(base)/$(bptype)/mixed/stats"

opt_dim:
	python opt_dim.py $(bptype) mixed
	pandoc -V geometry:landscape -o $(bptype)/mixed/opt_dim.pdf $(bptype)/mixed/opt_dim.txt
	$(DBU) upload $(bptype)/mixed/opt_dim.pdf "$(base)/$(bptype)/mixed"

growth_patterns_nocontrol:
	python growth_patterns_nocontrol.py $(bptype) 1 3 nodadbp
	python growth_patterns_nocontrol.py $(bptype) 1 3 dadbp
	python growth_patterns_nocontrol_plot.py $(bptype) 1 nodadbp
	python growth_patterns_nocontrol_plot.py $(bptype) 1 dadbp
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_nocontrolcbs_nodadbp_1_3.pdf $(bptype)/mixed/dim_1/growth_patterns_nocontrolcbs_nodadbp_1_3.txt
	pandoc -o $(bptype)/mixed/dim_1/growth_patterns_nocontrolcbs_dadbp_1_3.pdf $(bptype)/mixed/dim_1/growth_patterns_nocontrolcbs_dadbp_1_3.txt
	$(DBU) upload $(bptype)/mixed/dim_1/growth_patterns_nocontrolcbs_nodadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
	$(DBU) upload $(bptype)/mixed/dim_1/growth_patterns_nocontrolcbs_dadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
	$(DBU) upload $(bptype)/mixed/dim_1/growth_patterns_plot_nocontrolcbs_nodadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
	$(DBU) upload $(bptype)/mixed/dim_1/growth_patterns_plot_nocontrolcbs_dadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"

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
	$(DBU) upload $(bptype)/mixed/dim_1/growth_patterns_controlcbs_dadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
	$(DBU) upload $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_4.pdf "$(base)/$(bptype)/mixed/dim_1"
	$(DBU) upload $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
	$(DBU) upload $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_3_zm.pdf "$(base)/$(bptype)/mixed/dim_1"
	$(DBU) upload $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_1.pdf "$(base)/$(bptype)/mixed/dim_1"
	$(DBU) upload $(bptype)/mixed/dim_1/growth_patterns_controlcbs_nodadbp_1_3_zf.pdf "$(base)/$(bptype)/mixed/dim_1"
	$(DBU) upload $(bptype)/mixed/dim_1/growth_patterns_plot_controlcbs_nodadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
	$(DBU) upload $(bptype)/mixed/dim_1/growth_patterns_plot_controlcbs_dadbp_1_3.pdf "$(base)/$(bptype)/mixed/dim_1"
