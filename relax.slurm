#!/bin/bash                                                          
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=relax
#SBATCH --ntasks=32
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
# #SBATCH --mail-user=YOUR@EMAIL
# #SBATCH --mail-type=END

module purge
module load QuantumESPRESSO/7.3-foss-2023a

mpirun pw.x -inp relax.in > relax.out

