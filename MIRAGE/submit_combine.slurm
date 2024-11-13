#!/bin/bash -l
#SBATCH -t 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=dimml002@umn.edu 
#SBATCH --job-name combine
#SBATCH -o ./jobs/%j.combine.out
#SBATCH -e ./jobs/%j.combine.err
#SBATCH -p agsmall,aglarge

set -e

cd $SLURM_SUBMIT_DIR

micromamba activate impute

snakemake -s combine.smk \
    --profile slurm.PaintedPatches/ \
    --rerun-triggers mtime \
    --keep-going \
    --use-conda

