#!/bin/bash

# Loop through all slurm-*.out and scancel the corresponding jobs
for file in slurm-*.out
do
    jobid=$(echo $file | cut -d'-' -f2 | cut -d'.' -f1)
    echo "scancel $jobid"
    scancel $jobid
    #rm $file
done
