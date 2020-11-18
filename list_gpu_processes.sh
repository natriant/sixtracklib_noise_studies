#! /bin/bash

printf "\n"
for gpu_index in 0 1 2 3
do
echo "GPU-ID: $gpu_index (`nvidia-smi -i 0 --query-gpu=name --format=csv,noheader`)"
join -1 1 -2 3 \
<(nvidia-smi -i $gpu_index --query-compute-apps=pid,used_memory \
--format=csv \
| sed "s/ //g" | sed "s/,/ /g" \
| awk 'NR<=1 {print toupper($0)} NR>1 {print $0}' \
| sed "/\[NotSupported\]/d" \
| awk 'NR<=1{print $0;next}{print $0| "sort -k1"}') \
<(ps -a -o user,pgrp,pid,pcpu,pmem,time,command \
| awk 'NR<=1{print $0;next}{print $0| "sort -k3"}') \
| column -t
printf "\n"
done
