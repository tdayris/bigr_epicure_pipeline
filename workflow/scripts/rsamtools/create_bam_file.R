# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script create a BamFile object for future use

base::library(package = "Rsamtools", character.only = TRUE)
base::message("Library loaded")

bam_file <- Rsamtools::BamFile(
    file = base::as.character(x = snakemake@input[["bam"]]),
    index = base::as.character(x = snakemake@input[["bai"]])
)
base::message("Bam file created")

base::saveRDS(
    object = bam_file,
    file = base::as.character(x = snakemake@output["rds"])
)
base::message("Process over")

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()