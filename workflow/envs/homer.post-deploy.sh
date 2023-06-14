set -eo pipefail

perl "${CONDA_PREFIX}share/homer/configureHomer.pl" -install hg38 -keepScript
perl "${CONDA_PREFIX}share/homer/configureHomer.pl" -install mm10 -keepScript
