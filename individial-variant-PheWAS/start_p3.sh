#!/bin/bash

dir_name='stable'

log_name="p3-$(date '+%Y-%m-%d_%H-%M-%S').txt"

R --slave --no-restore \
--file=Step2_p3.r \
--args \
--jobs_path /pathto/jobs_num.Rdata \
--pheno_path /pathto/disease/ \
--step1_path /pathto/STEP_1/ \
--out_path /pathto/STEP_2/${dir_name}/ \
--step2.1_path /pathto/STEP_2/${dir_name}/step2.1/ \
--step2.2_path /pathto/STEP_2/${dir_name}/step2.2/ \
--log_path /pathto/STEP_2/log/${dir_name} \
--arrayid 25 \
--begin 141 \
--end 150 \
--stage train > "$log_name" 2>&1
