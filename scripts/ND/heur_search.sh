#!/bin/sh

# SLURM options:
#SBATCH --partition=htc
#SBATCH --array=0-4
#SBATCH --job=RH
#SBATCH --licenses=sps
#SBATCH --time=0-16
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/ND/heur_%A_%a.log

PROJECT=/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate
SCRIPTDIR="/sps/nemo/scratch/mpetro/Projects/PhD/SNSensitivityEstimate/scripts/ND"

echo "chose signal process: bb0nu_foil_bulk, bb0nuM1_foil_bulk, bb0nuM2_foil_bulk, RH037_foil_bulk"
SIGNAL="RH050_foil_bulk"
echo "sending job for signal: $SIGNAL"

# change signal parameter in SCRIPT via sed and copy to new file
sed -e "s|%SIGNAL|$SIGNAL|" \
     $SCRIPTDIR/heuristic.jl > $SCRIPTDIR/heuristic_$SIGNAL.jl

cd $PROJECT

module load julia

julia --project=$PROJECT -e "include(\"$SCRIPTDIR/heuristic_$SIGNAL.jl\")"
