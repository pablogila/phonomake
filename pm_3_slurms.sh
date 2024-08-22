#!/bin/bash

TEMPLATE="scf.slurm"
INPUT_FILE_KEYWORD="INPUT_FILE"
OUTPUT_FILE_KEYWORD="OUTPUT_FILE"
JOB_NAME_KEYWORD="JOB_NAME"
TRASH_FOLDER="slurms"

# Create a trash folder to later move the files if it doesn't exist
mkdir -p $TRASH_FOLDER
# Loop through all supercell-*.in files
for file in supercell-*.in
do
    # Generate unique job script name
    job_script="${file%.in}.slurm"
    # Replace placeholders in the template with actual file names
    sed "s|$INPUT_FILE_KEYWORD|$file|g; s|$OUTPUT_FILE_KEYWORD|${file%.in}.out|g; s|$JOB_NAME_KEYWORD|${file%.in}|g" $TEMPLATE > "$job_script"
    # Submit the job
    #source "$job_script"
    sbatch "$job_script"
    # Move the job script to the trash folder
    mv "$job_script" "$TRASH_FOLDER/"
done

