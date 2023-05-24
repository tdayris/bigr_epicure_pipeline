# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script calculates the local densities of a defined 
# sequence pattern (e.g. CpGs) 

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")


# Load libraries
base::library(package = "MEDIPS", character.only = TRUE)
base::library(package = "BSgenome.Hsapiens.UCSC.hg38", character.only = TRUE)
base::library(package = "BSgenome.Mmusculus.UCSC.mm10", character.only = TRUE)
base::message("Libraries loaded")

medips_sets <- base::readRDS(
    file = base::as.character(x = snakemake@input[["rds"]])
)
base::message("MEDIPS-sets loaded")

coupling <- MEDIPS::MEDIPS.couplingVector(
    pattern = base::as.character(x = snakemake@params[["pattern"]]),
    refObj = medips_sets[[1]]
)
base::message("MEDIPS-sets CpG density computed")


saveRDS(
    object = coupling,
    file = snakemake@output[["rds"]]
)
base::message("Process over")

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
sink(type = "message")
sink()