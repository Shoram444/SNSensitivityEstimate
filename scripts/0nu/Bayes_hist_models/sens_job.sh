#!/bin/sh

# SLURM options:
#SBATCH --partition=htc
#SBATCH --job=0nuM1_bayes
#SBATCH --licenses=sps
#SBATCH --array=0-9
#SBATCH --time=1-6
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --output=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/slurm_jobs/bayes_logs/bayes_%job_%A_%a.log


PROJECT=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate
SCRIPTDIR="/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/0nu/Bayes_hist_models/"

# bb0nu_foil_bulk, bb0nuM1_foil_bulk, bb0nuM2_foil_bulk, RH037_foil_bulk

SIGNAL="bb0nuM2_foil_bulk"
echo "signal process: $SIGNAL"
BINLOW=300
echo "BINLOW: $BINLOW"
BINHIGH=3500
echo "BINHIGH: $BINHIGH"
BINWIDTH=50
echo "BINWIDTH: $BINWIDTH"
MODE=sumE
echo "MODE: $MODE"
PRIOR=1e-2
echo "PRIOR: $PRIOR"

echo "sending job for signal: $SIGNAL"

# change signal parameter in SCRIPT via sed and copy to new file
sed -e "s|%SIGNAL|$SIGNAL|" \
    -e "s|%BINLOW|$BINLOW|" \
    -e "s|%BINHIGH|$BINHIGH|" \
    -e "s|%BINWIDTH|$BINWIDTH|" \
    -e "s|%MODE|$MODE|" \
    -e "s|%PRIOR|$PRIOR|" \
     $SCRIPTDIR/scratchpad.jl > $SCRIPTDIR/scripts/bayes_$SIGNAL.jl

cd $PROJECT

SCRIPT="$SCRIPTDIR/scripts/bayes_$SIGNAL.jl"

module load julia
export JULIA_NUM_THREADS=4

julia --threads=4 --project=$PROJECT -e "include(\"$SCRIPT\")"
# julia --threads=4 --project=$PROJECT -e "using Pkg; Pkg.instantiate(); include(\"$SCRIPT\")"


