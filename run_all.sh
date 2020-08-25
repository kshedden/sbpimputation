#!/bin/bash

if [ -z "$1" ]
then
    echo "Required arguments: [sbp|dbp]"
    exit 1
fi

python run_sbp_pca.py $1 1 controlcbs dadbp >& t1
python run_sbp_pca.py $1 1 nocontrolcbs dadbp >& t2
python run_sbp_pca.py $1 1 controlcbs nodadbp >& t3
python run_sbp_pca.py $1 1 nocontrolcbs nodadbp >& t4

python run_sbp_pca_sex.py $1 1 controlcbs dadbp female >& t5
python run_sbp_pca_sex.py $1 1 nocontrolcbs dadbp female >& t6
python run_sbp_pca_sex.py $1 1 controlcbs nodadbp female >& t7
python run_sbp_pca_sex.py $1 1 nocontrolcbs nodadbp female >& t8

python run_sbp_pca_sex.py $1 1 controlcbs dadbp male >& t9
python run_sbp_pca_sex.py $1 1 nocontrolcbs dadbp male >& t10
python run_sbp_pca_sex.py $1 1 controlcbs nodadbp male >& t11
python run_sbp_pca_sex.py $1 1 nocontrolcbs nodadbp male >& t12

python run_sbp_pca.py $1 2 controlcbs dadbp >& t13
python run_sbp_pca.py $1 2 nocontrolcbs dadbp >& t14
python run_sbp_pca.py $1 2 controlcbs nodadbp >& t15
python run_sbp_pca.py $1 2 nocontrolcbs nodadbp >& t16

