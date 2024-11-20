#!/bin/bash
#SBATCH --time=0-00:15:00
#SBATCH --partition=short
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=outfile_summary
module purge
module load R

### Run simulation code and benchmark the start/end times
start_time=$(date +%s)
Rscript summary_plots.R
finish_time=$(date +%s)

### Calculate and output total runtime
elapsed_time=$((finish_time  - start_time))
((sec=elapsed_time%60, elapsed_time/=60, min=elapsed_time%60, hrs=elapsed_time/60))
timestamp=$(printf "Total time taken - %d hours, %d minutes, and %d seconds." $hrs $min $sec)
echo $timestamp

