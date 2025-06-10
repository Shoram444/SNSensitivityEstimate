#!/bin/sh

# SLURM options:
#SBATCH --partition=flash
# SBATCH --array=0-0
#SBATCH --job=0nu
#SBATCH --licenses=sps
#SBATCH --time=00:20
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/slurm_jobs/nDim_logs/0nu_1_nDim_%A_%a.log

PROJECT=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate
SCRIPT="/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/ND/scratchpad.jl"


cd $PROJECT

module load julia

julia --project=$PROJECT -e "include(\"$SCRIPT\")"
