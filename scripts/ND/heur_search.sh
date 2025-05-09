#!/bin/sh

# SLURM options:
#SBATCH --partition=htc
#SBATCH --array=0-0
#SBATCH --job=heur
#SBATCH --licenses=sps
#SBATCH --time=1-0
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/ND/heur_%A_%a.log

PROJECT=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate
SCRIPT="/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/ND/heuristic.jl"


cd $PROJECT

module load julia

julia --project=$PROJECT -e "include(\"$SCRIPT\")"
