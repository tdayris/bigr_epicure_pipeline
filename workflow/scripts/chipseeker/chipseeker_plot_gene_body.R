# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script gene body coverage using CHiPSeeker

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")
base::message("Logging defined")

# Load libraries
base::library(package = "ChIPseeker", character.only = TRUE)
base::library(
    package = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    character.only = TRUE
)
base::library(
    package = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    character.only = TRUE
)
base::message("Libraries loaded")

ranges <- NULL
if ("bed" %in% base::names(x = snakemake@input)) {
    bed_path <- base::as.character(x = snakemake@input[["bed"]])
    ranges <- ChIPseeker::readPeakFile(
        peakfile = bed_path, as = "GRanges", header = FALSE
    )
} else if ("ranges" %in% base::names(x = snakemake@input)) {
    ranges <- base::readRDS(
        file = base::as.character(x = snakemake@input[["ranges"]])
    )
}
tag_matrix <- base::readRDS(
    file = base::as.character(x = snakemake@input[["tagmatrix"]])
)
base::message("File list acquired")


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
if ("organism" %in% base::names(x = snakemake@params)) {
    if (snakemake@params[["organism"]] == "mm10") {
        organism <- TxDb.Mmusculus.UCSC.mm10.knownGene
    }
}

# Build plot
png(
  filename = snakemake@output[["png"]],
  width = 1024,
  height = 768,
  units = "px",
  type = "cairo"
)

ChIPseeker::plotPeakProf2(
    peak = ranges,
    tagMatrix = tag_matrix,
    upstream = rel(0.2),
    downstream = rel(0.2),
    weightCol = "V5",
    ignore_strand = TRUE,
    conf = 0.95,
    by = "gene",
    type = "body",
    TxDb = txdb,
    facet = "row",
    nbin = 800
)

dev.off()
base::message("Plot saved")

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()