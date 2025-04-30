#!/bin/sh

# SLURM options:
#SBATCH --partition=htc
#SBATCH --job=bayes_M2
#SBATCH --licenses=sps
#SBATCH --array=0-19
#SBATCH --time=0-40
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --output=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/slurm_jobs/bayes_logs/bayes_%A_%a.log


PROJECT=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate
SCRIPT="/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/0nu/Bayes_hist_models/scratchpad.jl"

cd $PROJECT

module load julia
export JULIA_NUM_THREADS=4

julia --threads=4 --project=$PROJECT -e "include(\"$SCRIPT\")"
# julia --threads=4 --project=$PROJECT -e "using Pkg; Pkg.instantiate(); include(\"$SCRIPT\")"