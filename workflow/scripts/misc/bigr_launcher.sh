#!/usr/bin/bash
set -euiop pipefail
shopt -s nullglob

#SBATCH --job-name="BiGR_EpiCure_Pipeline"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --partition=mediumq

########################
## EDIT THE FOLLOWING ##
########################

# Path to your working directory
WORK_DIR=""

# Sequencing protocol type: "atac-seq", "chip-seq", "medip-seq", "og-seq"
PROTOCOL="chip-seq"

# Keep the pipeline version homogen
# among multiple sequencing in a same
# research project
PIPELINE_VERSION="0.14.1"

# Organism name: "homo_sapiens" or "mus_musculus"
ORGANISM="homo_sapiens"

# Library type: single or paired
LIBRARY="paired"

########################################
##        DO NOT EDIT BELOW           ##
## Unless you know what you are doing ##
########################################

echo "Working in: ${WORK_DIR}"
WORK_DIR=$(readlink -e "${WORK_DIR:?}")
cd "${WORK_DIR}" || exit 1

# Activating shell environment
echo "Activating conda..."
source /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/bigr_snakemake/etc/profile.d/conda.sh
conda activate "/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/shared_install/bigr_epicure_pipeline"

# Exporting required IO variables
echo "Exporting environment variables..."
export SNAKEMAKE_OUTPUT_CACHE="/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/snakemake_cache/"

# Building IO architecture
echo "Building repository architecture if missing..."
mkdir --parent --verbose logs/slurm tmp data_{output,input} reference/{blacklist,xenome/index}
# cd tmp || exit 1

# Deploy workflow
if [ ! -f "workflow/Snakefile" ]; then
    echo "Deploying Snakemake-workflow..."
    snakedeploy deploy-workflow "https://github.com/tdayris/bigr_epicure_pipeline" . --tag "${PIPELINE_VERSION}"
    mv "config/design.tsv" "config/example.design.tsv"
else
    echo "Updating Snakemake-workflow verison if needed..."
    sed "s|tag=\"[^\"]\"|tag=\"${PIPELINE_VERSION}\"|g" "workflow/Snakefile"
fi

if [ ! -f "workflow/bigr_launcher.sh" ]; then
    echo "Saving local copy of this launcher script..."
    LAUNCHER_PATH=$(readlink -e "$0")
    cp --verbose "${LAUNCHER_PATH}" "workflow/bigr_launcher.sh"
fi
    

# Edit configuration file if needed
if [ "${ORGANISM}" == "mus_musculus" ]; then
    echo "Editing configuration in order to match Mouse datasets"
    sed --in-place 's|organism: homo_sapiens|organism: mus_musculus|g' "config/config.yaml"
    sed --in-place 's|build: GRCh38|build: GRCm38|g' "config/config.yaml"
    sed --in-place 's|release: 109|release: 103|g' "config/config.yaml"
fi

# Update sequencing protocol
echo "Updating sequencing protocol if needed..."
sed --in-place \
    "s|protocol: atac-seq|protocol: ${PROTOCOL}|g" \
    "config/config.yaml"

# Use BiGR Annotations
echo "Linking available annotations (if any) ..."
sed --in-place \
    's|fastq_screen_config: config/fastq_screen.conf|fastq_screen_config: config/fastq_screen_bigr.conf|g' \
    "config/config.yaml"

if [ "${ORGANISM}" == "homo_sapiens" ]; then
    ln --symbolic --force --relative --verbose \
        "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta" \
        "reference/homo_sapiens.GRCh38.109.fasta"
  
    ln --symbolic --force --relative --verbose \
        "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta.fai" \
        "reference/homo_sapiens.GRCh38.109.fasta.fai"

    ln --symbolic --force --relative --verbose \
        "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.dict" \
        "reference/homo_sapiens.GRCh38.109.dict"
 
    ln --symbolic --force --relative --verbose \
        "/mnt/beegfs/database/bioinfo/Index_DB/GTF/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.gtf" \
        "reference/homo_sapiens.GRCh38.109.gtf"

    ln --symbolic --force --relative --verbose \
        "/mnt/beegfs/database/bioinfo/Index_DB/Bowtie/2.5.1/GRCh38.109/bowtie2_index/" \
        "reference/"

    ln --symbolic --force --relative --verbose \
        "/mnt/beegfs/database/bioinfo/Index_DB/2bit/GRCh38.109/homo_sapiens.GRCh38.109.2bit" \
        "reference/homo_sapiens.GRCh38.109.2bit"

    ln --symbolic --force --relative --verbose \
        "/mnt/beegfs/database/bioinfo/Index_DB/blacklist/hg38.blacklist.merged.bed" \
        "reference/blacklist/homo_sapiens.GRCh38.merged.109.bed.gz"

    find "/mnt/beegfs/database/bioinfo/Index_DB/xenome/GRCh38_GRCm38/Ensembl/r99/DNA/" -type f | while read FILE; do
        FILE_NAME=$(basename "${FILE/GRCh38_GRCm38_r99_DNA/pdx-both}")
        ln --symbolic --force --relative --verbose \
            "${FILE}" \
            "${FILE_NAME}"
    done

elif [ "${ORGANISM}" == "mus_musculus" ]; then
    echo "Mouse references not yet available"
fi

# Build design file if missing
if [ ! -f "config/design.tsv" ]; then
    echo "Building design file..."
    if [ "${LIBRARY}" == "paired" ]; then
        echo "Looking for fastq files..."
        find "${WORK_DIR}/data_input/" -type f -name "*q.gz" | \
            sort | \
            uniq | \
            paste - - | \
            while read R1 R2; do echo -e "$(basename ${R1})\t${R1}\t${R2}"; done | \
            awk 'BEGIN{FS="\t"; print "Sample_id" FS "Upstream_file" FS "Downstream_file"} {print $0}' \
        > "config/design.tsv"
    else
        echo "Looking for fastq files..."
        find "${WORK_DIR}/data_input/" -type f -name "*q.gz" | \
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

# Launching pipeline
echo "Launching pipeline..."
snakemake --snakefile 'workflow/Snakefile' \
          --profile '/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/profiles/slurm-web/' \
          --cache \
          --rerun-incomplete \
          --nt "$@"