#!/bin/sh

# SLURM options:
#SBATCH --partition=htc
#SBATCH --licenses=sps
#SBATCH --array=0-9
#SBATCH --time=0-42
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/slurm_jobs/bayes_%A.log

JULIA=~/.julia/juliaup/julia-1.10.0+0.x64.linux.gnu/bin
PROJECT=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate
SCRIPT="/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/0nu/Bayes_hist_models/scratchpad.jl"

$JULIA/julia --project=$PROJECT -e "using Pkg; Pkg.instantiate(); include(\"$SCRIPT\")"