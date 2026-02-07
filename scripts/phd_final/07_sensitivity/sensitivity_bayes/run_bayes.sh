#!/bin/sh

# SLURM options:
#SBATCH --partition=htc
#SBATCH --job=0nu2
#SBATCH --licenses=sps
#SBATCH --array=0-9
#SBATCH --time=0-2 #1-6
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --output=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/phd_final/07_sensitivity/sensitivity_bayes/logs/bayes_%job_%A_%a.log


PROJECT=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate
SCRIPTDIR="/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/phd_final/07_sensitivity/sensitivity_bayes/"

# bb0nu_foil_bulk, bb0nuM1_foil_bulk, bb0nuM2_foil_bulk, RH037_foil_bulk

SIGNAL="bb0nu_foil_bulk"
echo "signal process: $SIGNAL"
BINLOW=300
echo "BINLOW: $BINLOW"
BINHIGH=3500
echo "BINHIGH: $BINHIGH"
BINWIDTH=100
echo "BINWIDTH: $BINWIDTH"
MODE=sumE
echo "MODE: $MODE"
PRIOR=1e-3
echo "PRIOR: $PRIOR"
RADON_TAG=2
echo "RADON_TAG: $RADON_TAG"

echo "sending job for signal: $SIGNAL"

# change signal parameter in SCRIPT via sed and copy to new file
sed -e "s|%SIGNAL|$SIGNAL|" \
    -e "s|%BINLOW|$BINLOW|" \
    -e "s|%BINHIGH|$BINHIGH|" \
    -e "s|%BINWIDTH|$BINWIDTH|" \
    -e "s|%MODE|$MODE|" \
    -e "s|%PRIOR|$PRIOR|" \
    -e "s|%RADON_TAG|$RADON_TAG|" \
     $SCRIPTDIR/bayes_run_slurm.jl > $SCRIPTDIR/scripts/bayes_${SIGNAL}_radon${RADON_TAG}.jl

cd $PROJECT

SCRIPT="$SCRIPTDIR/scripts/bayes_${SIGNAL}_radon${RADON_TAG}.jl"
module load julia
export JULIA_NUM_THREADS=4

julia --threads=4 --project=$PROJECT -e "include(\"$SCRIPT\")"
# julia --threads=4 --project=$PROJECT -e "using Pkg; Pkg.instantiate(); include(\"$SCRIPT\")"


