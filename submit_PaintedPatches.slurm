#!/bin/bash -l
#SBATCH -t 6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=dimml002@umn.edu 
#SBATCH --job-name PaintedPatches
#SBATCH -o ./jobs/%j.PaintedPatches.out
#SBATCH -e ./jobs/%j.PaintedPatches.err
#SBATCH -p amdsmall,amdlarge

set -e

cd $SLURM_SUBMIT_DIR

micromamba activate impute

snakemake -s PaintedPatches.smk \
    --profile slurm.PaintedPatches/ \
    --configfile config.yaml \
    --rerun-triggers mtime \
    --keep-going \
    --use-conda

