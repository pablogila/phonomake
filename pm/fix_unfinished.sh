#!/bin/bash

# This bash script checks which supercell-* calculations did not finish, and re-runs them
# If the file stopped due to the time limit, it edits the supercell-*.in file to set restart_mode to restart
# If it stopped abruptly, it removes the old output file and starts again
# If the tail of the file does not contain "JOB DONE.", it re-runs the file, with the new memory limit
# The scf.slurm file is used as a template for the new slurm files

# Prompt the user for the new memory value
echo "Let's fix the unfinished 'supercell-*' calculations!"
read -p "New memory limit (int, GB, 0 to keep old values):  " MEMORY

SLURM_TEMPLATE="scf.in"
JOB_NAME_KEYWORD="JOB_NAME"
INPUT_FILE_KEYWORD="INPUT_FILE"
OUTPUT_FILE_KEYWORD="OUTPUT_FILE"
TIME_EXCEEDED="Maximum CPU time exceeded"
RESTART_MODE_KEYWORD="restart_mode"
RESTART_MODE="restart_mode = 'restart'"
START_MODE="restart_mode = 'from_scratch'"
JOB_DONE="JOB DONE."
SLURMS_FOLDER="slurms"

for file in supercell-*.out; do

    # Check that the file exists
    if [ ! -f "$file" ]; then
        echo "No $file file was found. Wtf are you doing."
        return
    fi

    tail=$(tail -n 200 "$file")
    # Skip the file if it ended properly
    if [[ $tail == *"$JOB_DONE"* ]]; then
        if [[ $tail != *"$TIME_EXCEEDED"* ]]; then
            echo ">>> Skipped finished $file"
            continue
        fi
    fi

    sbatch_file="${file%.out}.slurm"
    cp scf.slurm "$sbatch_file"
    sed -i "s/$JOB_NAME_KEYWORD/${file%.out}/" "$sbatch_file"
    sed -i "s/INPUT_FILE/${file%.out}.in/" "$sbatch_file"
    sed -i "s/OUTPUT_FILE/$file/" "$sbatch_file"

    # Increase the memory limit to the new value
    if [ $MEMORY -ne 0 ]; then
        sed -i "s/#SBATCH --mem=.*$/#SBATCH --mem=${MEMORY}G/" "$sbatch_file"
    fi

    # Restart the file if it stopped due to the time limit, else remove the old one
    if [[ $tail == *"$TIME_EXCEEDED"* ]]; then
        sed -i "s/$RESTART_MODE_KEYWORD.*$/$RESTART_MODE/" "${file%.out}.in"
        echo "··· Resumed unfinished $file"
    else
        sed -i "s/$RESTART_MODE_KEYWORD.*$/$START_MODE/" "${file%.out}.in"
        rm -f $file
        echo "--- Removed faulty $file"
    fi

    sbatch "$sbatch_file"
    #source "$sbatch_file"
    mkdir -p $SLURMS_FOLDER
    mv "$sbatch_file" $SLURMS_FOLDER
done
