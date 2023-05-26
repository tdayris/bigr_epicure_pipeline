#!/usr/bin/bash
#SBATCH --job-name="BiGR_EpiCure_Pipeline"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --partition=mediumq


conda shell.bash activate /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/shared_install/bigr_epicure_pipeline

export SNAKEMAKE_OUTPUT_CACHE="/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/snakemake_cache/"

mkdir --parent --verbose logs/slurm

snakemake --snakefile 'workflow/Snakefile' \
          --cache \
          --profile '/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.2.0/bigr_pipelines/common/profiles/slurm'
