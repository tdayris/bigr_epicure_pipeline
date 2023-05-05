set -e
set -u
set -o pipefail
set -p

if [[ $(conda info --envs | grep -q bigr_epicure_pipeline) -eq 1 ]]; then 
    mamba env create -f ../workflow/envs/bigr_epicure_pipeline.yaml
fi

#Overload general cache in order to test pipeline deployment
TMP="temp"
SNAKEMAKE_OUTPUT_CACHE="./cache"
export TMP SNAKEMAKE_OUTPUT_CACHE

# Run resources donload/indexing
conda shell.bash activate bigr_epicure_pipeline
snakemake -s ../workflow/Snakefile -c 1 --use-conda --cache --configfile ../config/config.install.yaml -n

# Delete output after testing
# snakemake -s ../workflow/Snakefile -c 1 --use-conda --cache --configfile ../config/config.install.yaml --delete-all-output
