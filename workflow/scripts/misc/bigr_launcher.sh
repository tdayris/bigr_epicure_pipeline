#!/usr/bin/bash
#SBATCH --job-name="BiGR_EpiCure_Pipeline"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --partition=mediumq


conda shell.bash activate bigr_epicure_pipeline

mkdir --parent --verbose logs/slurm

snakemake --snakefile 'workflow/Snakefile' \
          --jobs 30 \
          --local-cores 2 \
          --restart-times 1 \
          --max-jobs-per-second 1 \
          --max-status-checks-per-second 1 \
          --configfile 'config/config.yaml' \
          --keep-going \
          --reason \
          --printshellcmds \
          --rerun-triggers 'mtime' \
          --cluster-cancel 'scancel' \
          --cluster 'scripts/jobscripts/bigr_slurm_submit.py' \
          --jobscript 'scripts/jobscripts/jobscript.sh' \
          --cluster-status 'scripts/jobscripts/bigr_slurm_status.py' \
          --jobname "{name}.{jobid}.snakejob.sh" \
          --use-conda \
          --conda-prefix '/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/shared_install/' \
          --shadow-prefix 'tmp'
