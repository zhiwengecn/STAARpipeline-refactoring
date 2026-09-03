#!/bin/bash

prefix="stable"
log_name="p1-$(date '+%Y-%m-%d_%H-%M-%S').txt"

R --slave --no-restore \
--file=/pathto/${prefix}/2_Single_p1_gzw.r \
--args \
--adgs /pathto/agds_dir.Rdata \
--null /pathto/obj_nullmodel.Rdata \
--jobs /pathto/jobs_num.Rdata \
--p1_path /pathto/STEP_2/P1/${prefix}/ \
--log_path /pathto/logs/${prefix}/ \
--array_start 2 \
--array_end 2 \
--subset_num 5000 > "$log_name" 2>&1
