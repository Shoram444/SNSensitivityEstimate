#!/bin/sh

# SLURM options:
#SBATCH --partition=htc
#SBATCH --job=rh020
#SBATCH --licenses=sps
#SBATCH --array=0-0
#SBATCH --time=1-6
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --output=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/phd_final/07_sensitivity/sensitivity_nd/nd_cc_runs/logs/nd_%job_%A_%a.log


PROJECT=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate
SCRIPTDIR="/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/phd_final/07_sensitivity/sensitivity_nd/nd_cc_runs/"

# bb0nu_foil_bulk, bb0nuM1_foil_bulk, bb0nuM2_foil_bulk, RH050_foil_bulk, _foil_bulk, Nnubb1500keV_foil_bulk

SIGNAL="Nnubb500keV_foil_bulk"
echo "signal process: $SIGNAL"
HOURS=22
echo "HOURS: $HOURS"
RADON_TAG=1
echo "RADON_TAG: $RADON_TAG"
SIDE="sameSide" # "both" or "sameSide", "oppositeSide" 

echo "sending job for signal: $SIGNAL"

# change signal parameter in SCRIPT via sed and copy to new file
sed -e "s|%SIGNAL|$SIGNAL|" \
    -e "s|%HOURS|$HOURS|" \
    -e "s|%RADON_TAG|$RADON_TAG|" \
     $SCRIPTDIR/nd_script.jl > $SCRIPTDIR/scripts/nd_${SIGNAL}_radon${RADON_TAG}.jl

cd $PROJECT

SCRIPT="$SCRIPTDIR/scripts/nd_${SIGNAL}_radon${RADON_TAG}.jl"
module load julia
export JULIA_NUM_THREADS=8

julia --threads=8 --project=$PROJECT -e "include(\"$SCRIPT\")"
# julia --threads=4 --project=$PROJECT -e "using Pkg; Pkg.instantiate(); include(\"$SCRIPT\")"


