# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script filters counted reads with csaw

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")


# Load libraries
base::library(package = "ChIPseeker", character.only = TRUE)
base::library(package = "org.Hs.eg.db", character.only = TRUE)
base::library(package = "org.Mm.eg.db", character.only = TRUE)
base::message("Libraries loaded")