#!/bin/sh

# SLURM options:
#SBATCH --partition=htc
#SBATCH --job=rh_oppo
#SBATCH --licenses=sps
#SBATCH --array=0-0
#SBATCH --time=1-0
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --output=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/phd_final/07_sensitivity/sensitivity_nd/nd_cc_runs/logs/nd_%job_%A_%a.log


PROJECT=/Users/maros.petro/Work/Phd/SNSensitivityEstimate #/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate
SCRIPTDIR=/Users/maros.petro/Work/Phd/SNSensitivityEstimate/scripts/phd_final/07_sensitivity/sensitivity_nd/nd_cc_runs/ #"/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/phd_final/07_sensitivity/sensitivity_nd/nd_cc_runs/"

# bb0nu_foil_bulk, bb0nuM1_foil_bulk, bb0nuM2_foil_bulk, RH050_foil_bulk, _foil_bulk, Nnubb1500keV_foil_bulk

SIGNAL="bb0nu_foil_bulk"
echo "signal process: $SIGNAL"
HOURS=2
echo "HOURS: $HOURS"
RADON_TAG=2
echo "RADON_TAG: $RADON_TAG"
SIDE="sameSide" # "both" or "sameSide", "oppositeSide" 

echo "sending job for signal: $SIGNAL"

# change signal parameter in SCRIPT via sed and copy to new file
sed -e "s|%SIGNAL|$SIGNAL|" \
    -e "s|%HOURS|$HOURS|" \
    -e "s|%RADON_TAG|$RADON_TAG|" \
    -e "s|%SIDE|$SIDE|" \
     $SCRIPTDIR/nd_script.jl > $SCRIPTDIR/scripts/nd_${SIGNAL}_radon${RADON_TAG}_side${SIDE}.jl

cd $PROJECT

SCRIPT="$SCRIPTDIR/scripts/nd_${SIGNAL}_radon${RADON_TAG}_side${SIDE}.jl"
# module load julia
# export JULIA_NUM_THREADS=4

# julia --threads=4 --project=$PROJECT -e "include(\"$SCRIPT\")"

