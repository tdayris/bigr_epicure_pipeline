set -e

# bash pipeline_lint.sh

declare -a protocols
protocols=(
    "chip-seq"
    "atac-seq"
    "oxidip-seq"
    "medip-seq"
    "cutntag"
    "cutnrun"
)

for PROTOCOL in "${protocols[@]}"; do
    echo "Editing protocol as: ${PROTOCOL}"
    sed -i "s|protocol: .*|protocol: '${PROTOCOL}'|g" config/config.yaml

    echo "Building Rulegraph"
    snakemake -c 1 \
              --use-conda \
              -s ../workflow/Snakefile \
              --configfile ../config/config.yaml \
              --rulegraph | dot -T png > "../docs/rulegraphs/${PROTOCOL}.rulegraph.png"
done