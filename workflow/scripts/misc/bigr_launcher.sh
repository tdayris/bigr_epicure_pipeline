#!/usr/bin/bash
set -euiop pipefail
shopt -s nullglob

#SBATCH --job-name="BiGR_EpiCure_Pipeline"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --partition=mediumq

# Activating shell environment
echo "Activating conda..."
source /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/bigr_snakemake/etc/profile.d/conda.sh
conda activate "/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/shared_install/bigr_epicure_pipeline"

# Exporting required IO variables
echo "Exporting environment variables..."
export SNAKEMAKE_OUTPUT_CACHE="/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/snakemake_cache/"
export PIPELINE_VERSION="0.16.4"

# Building IO architecture
echo "Building repository architecture if missing..."
mkdir --parent --verbose logs/slurm tmp data_{output,input}

# Deploy workflow
if [ ! -f "workflow/Snakefile" ]; then
    echo "Deploying Snakemake-workflow..."
    snakedeploy deploy-workflow "https://github.com/tdayris/bigr_epicure_pipeline" . --tag "${PIPELINE_VERSION}"
    mv "config/design.tsv" "config/example.design.tsv"
fi

if [ ! -f "workflow/bigr_launcher.sh" ]; then
    echo "Saving local copy of this launcher script..."
    LAUNCHER_PATH=$(readlink -e "$0")
    cp --verbose "${LAUNCHER_PATH}" "workflow/bigr_launcher.sh"
fi

# Edit configuration file
if [ -f "config/config_BiGR_Flamingo.yaml" ]; then
    mv --verbose "config/config_BiGR_Flamingo.yaml" "config/config.yaml"
fi


# Build design file if missing
if [ ! -f "config/design.tsv" ]; then
    if [ -f "design.tsv" ]; then
        mv --verbose "design.tsv" "config/design.tsv"
    elif [ "${LIBRARY}" == "paired" ]; then
        echo "Looking for pairs of fastq files..."
        find "data_input/" -type f -name "*q.gz" | \
            sort | \
            uniq | \
            paste - - | \
            while read R1 R2; do echo -e "$(basename ${R1})\t${R1}\t${R2}"; done | \
            awk 'BEGIN{FS="\t"; print "Sample_id" FS "Upstream_file" FS "Downstream_file"} {print $0}' \
        > "config/design.tsv"
    else
        echo "Looking for fastq files..."
        find "data_input/" -type f -name "*q.gz" | \
            sort | \
            uniq | \
            while read R1; do echo -e "$(basename ${R1})\t${R1}"; done | \
            awk 'BEGIN{FS="\t"; print "Sample_id" FS "Upstream_file"} {print $0}' \
        > "config/design.tsv"
    fi
fi

# Dealing with windows/mac inputs
if [ -f "config/design.tsv" ]; then
    echo "Removing possible Windows complex end of lines..."
    sed -i 's|\r||g' config/design.tsv
fi

echo "Updating Snakemake-workflow verison if needed..."
sed -i "s|tag=\"[^\"]\"|tag=\"${PIPELINE_VERSION}\"|g" "workflow/Snakefile"

# Launching pipeline
echo "Launching pipeline..."
snakemake --snakefile 'workflow/Snakefile' \
            --profile '/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/profiles/slurm-web/' \
            --cache \
            --rerun-incomplete \
              --nt "$@"



