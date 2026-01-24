#!/bin/sh

# SLURM options:
#SBATCH --partition=htc
#SBATCH --array=0-4
#SBATCH --job=M2
#SBATCH --licenses=sps
#SBATCH --time=1-0
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/RH050_nDim_%A_%a.log

PROJECT=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate
SCRIPT="/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/phd_final/07_sensitivity/sensitivity_nd/nd_cc_runs/RH050_run.jl"

cd $PROJECT

module load julia

julia --project=$PROJECT -e "include(\"$SCRIPT\")"
