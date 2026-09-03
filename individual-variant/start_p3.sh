#!/bin/bash

prefix="stable"

R --slave --no-restore \
--file=/pathto/${prefix}/2_Single_p3_gzw.r \
--args \
--p2_path /pathto/STEP_2/P2/${prefix}/ \
--out_path /pathto/output/${prefix}/ \
--log_path /pathto/logs/${prefix}/ \
--array_id 2 \
--subset_num 5000 \
--chunk_size 100 \
--summary_file /pathto/STEP_2/P1/${prefix}/Analysis_Meta_Summary.txt
