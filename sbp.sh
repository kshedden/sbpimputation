#!/bin/bash

echo HT $1:HAZ $1:BAZ $1:WAZ $1:WT $1:BMI $1 | rush -D ":" -k 'python sbp_poly.py {}'

