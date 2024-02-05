#!/bin/bash -l
#SBATCH -t 12:00:00
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

micromamba activate combine_vcfs

snakemake -s combine.smk \
    --profile slurm.combine/ \
    --rerun-triggers mtime \
    --keep-going \
    --use-conda

