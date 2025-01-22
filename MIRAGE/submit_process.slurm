#!/bin/bash -l
#SBATCH -t 6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=dimml002@umn.edu 
#SBATCH --job-name process
#SBATCH -o ./jobs/%j.process.out
#SBATCH -e ./jobs/%j.process.err
#SBATCH -p amdsmall,amdlarge

set -e

cd $SLURM_SUBMIT_DIR

source /users/7/dimml002/miniforge3/etc/profile.d/conda.sh
conda activate impute

snakemake -s process.smk \
    --profile slurm.MIRAGE/ \
    --configfile process_config.yaml \
    --rerun-triggers mtime \
    --keep-going \
    --use-conda

