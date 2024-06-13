#!/bin/bash
#SBATCH --time=0-01:45:00
#SBATCH --mem-per-cpu=32GB
#SBATCH --output=outfile
module purge
module load R
Rscript generate_data.R
Rscript only_local.R

