set -eo pipefail

# Install organigme data
perl "${CONDA_PREFIX}/share/homer/configureHomer.pl" -install human-o -keepScript
perl "${CONDA_PREFIX}/share/homer/configureHomer.pl" -install mouse-o -keepScript

# Install genome data
perl "${CONDA_PREFIX}/share/homer/configureHomer.pl" -install hg38 -keepScript
perl "${CONDA_PREFIX}/share/homer/configureHomer.pl" -install mm10 -keepScript

# Install promoter data
perl "${CONDA_PREFIX}/share/homer/configureHomer.pl" -install human-p -keepScript
perl "${CONDA_PREFIX}/share/homer/configureHomer.pl" -install mouse-p -keepScript
