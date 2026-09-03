#!/bin/bash

dir_name='stable'

log_name="p1-$(date '+%Y-%m-%d_%H-%M-%S').txt"

R --slave --no-restore \
--file=Step2_p1.r \
--args \
--jobs_path /pathto/jobs_num.Rdata \
--adgs /pathto/agds_dir.Rdata \
--log_path /pathto/STEP_2/log/${dir_name} \
--pheno_path /pathto/disease/ \
--step1_path /pathto/STEP_1/ \
--out_path /pathto/STEP_2/${dir_name}/ \
--step2.1_path /pathto/STEP_2/${dir_name}/step2.1/ \
--begin 141 \
--end 150 \
--stage train \
--arrayid 25 > "$log_name" 2>&1
