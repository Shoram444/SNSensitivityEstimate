#!/bin/sh

# SLURM options:
#SBATCH --partition=htc
#SBATCH --job=bayes_RH
#SBATCH --licenses=sps
#SBATCH --array=0-19
#SBATCH --time=1-2
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/slurm_jobs/bayes_logs/bayes_%A_%a.log
#SBATCH --mail-user=maros.petro@cvut.cz
#SBATCH --mail-type=END


JULIA=~/.julia/juliaup/julia-1.10.0+0.x64.linux.gnu/bin
PROJECT=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate
SCRIPT="/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/0nu/Bayes_hist_models/scratchpad.jl"

$JULIA/julia --project=$PROJECT -e "using Pkg; Pkg.instantiate(); include(\"$SCRIPT\")"