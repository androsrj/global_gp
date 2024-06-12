#!/bin/bash
#SBATCH --time=0-02:00:00
#SBATCH --mem-per-cpu=32GB
#SBATCH --output=outfile
module purge
module load R
Rscript generate_data.R
Rscript full_gp.R

