#!/bin/sh

# SLURM options:
#SBATCH --partition=htc
#SBATCH --array=0-4
#SBATCH --job=M2
#SBATCH --licenses=sps
#SBATCH --time=2-0
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/M2_nDim_%A_%a.log

PROJECT=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate
SCRIPT="/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/phd_final/07_sensitivity/sensitivity_nd/nd_cc_runs/nu0bbM2_run.jl"


cd $PROJECT

module load julia

julia --project=$PROJECT -e "include(\"$SCRIPT\")"
