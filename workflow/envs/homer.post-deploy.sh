set -eo pipefail

# Install genome data
perl "${CONDA_PREFIX}share/homer/configureHomer.pl" -install hg38 -keepScript
perl "${CONDA_PREFIX}share/homer/configureHomer.pl" -install mm10 -keepScript

# Install promoter data
perl "${CONDA_PREFIX}share/homer/configureHomer.pl" -install human-p -keepScript
perl "${CONDA_PREFIX}share/homer/configureHomer.pl" -install mouse-p -keepScript