#!/bin/bash                                                          
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=JOB_NAME
#SBATCH --ntasks=64
#SBATCH --time=06:00:00
#SBATCH --mem=100G
# #SBATCH --mail-user=YOUR@EMAIL
# #SBATCH --mail-type=END

module purge
module load QuantumESPRESSO/7.2-foss-2023a

#echo "pw.x -inp INPUT_FILE > OUTPUT_FILE"
mpirun pw.x -inp INPUT_FILE > OUTPUT_FILE

