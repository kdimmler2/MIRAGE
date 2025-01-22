#!/bin/bash -l
#SBATCH -t 6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dimml002@umn.edu
#SBATCH --job-name missing
#SBATCH -o %j.missing.out
#SBATCH -e %j.missing.err
#SBATCH -p msismall,msilarge

source /users/7/dimml002/miniforge3/etc/profile.d/conda.sh
conda activate general_tools

bcftools filter -Oz -o no_missing.IDs.no_multis.pass.joint_call.goldenPath.20241024.vcf.gz -e 'F_MISSING > 0.0' /scratch.global/dimml002/process_vcf/results/split_multis/IDs.no_multis.pass.joint_call.goldenPath.20241024.vcf.gz
