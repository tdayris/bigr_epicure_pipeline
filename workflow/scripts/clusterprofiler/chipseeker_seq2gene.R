# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script converts peaks to a ClusterProfiler objects using CHiPSeeker

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")
base::message("Logging defined")

# Load libraries
base::library(package = "ChIPseeker", character.only = TRUE)
base::library(package = "TxDb.Hsapiens.UCSC.hg38.knownGene", character.only = TRUE)
base::library(package = "TxDb.Mmusculus.UCSC.mm10.knownGene", character.only = TRUE)
base::message("Libraries loaded")

peaks <- base::readRDS(
    file = base::as.character(x = snakemake@input[["peaks"]])
)
base::message("Peaks loaded")


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
if ("organism" %in% base::names(x = snakemake@params)) {
    if (snakemake@params[["organism"]] == "mm10") {
        organism <- TxDb.Mmusculus.UCSC.mm10.knownGene
    }
}


# Build seq2gene object
genes <- ChIPseeker::seq2gene(
    seq = peaks, 
    tssRegion = c(-1000, 1000),   
    flankDistance = 3000, 
    TxDb = txdb, 
    sameStrand = FALSE
)

# Save object
base::saveRDS(
    file = base::as.character(x = snakemake@output[["rds"]]),
    object = genes
)

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()