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

ncpus <- base::as.numeric(snakemake@threads[[1]])
mc.cores <- ncpus
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
base::message("Peaks/Ranges loaded")
tag_matrix <- base::readRDS(
    file = base::as.character(x = snakemake@input[["tagmatrix"]])
)
upstream <- base::attr(x = tag_matrix, which = "upsteam")
downstream <- base::attr(x = tag_matrix, which = "downstream")
base::message("Tagmatrix loaded")

# The default resample value
resample <- 500
window_number <- base::length(tag_matrix)

# Search all prime factors for the given number of windows
prime_factors <- function(x) {
  factors <- c()
  last_prime <- 2
  while(x >= last_prime){
    if (! x %% last_prime) {
      factors <- c(factors, last_prime)
      x <- x / last_prime
      last_prime <- last_prime - 1
      }
    last_prime <- last_prime + 1
  }
  base::return(factors)
}
primes <- prime_factors(x = window_number)
# Closest prime to defautl resample value
resample <- primes[which(abs(primes - resample) == min(abs(primes - resample)))]
base::message("Resampling value: ", resample)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
if ("organism" %in% base::names(x = snakemake@params)) {
    if (snakemake@params[["organism"]] == "mm10") {
        organism <- TxDb.Mmusculus.UCSC.mm10.knownGene
    }
}
base::message("Organism library available")

# Build plot
png(
  filename = snakemake@output[["png"]],
  width = 1024,
  height = 768,
  units = "px",
  type = "cairo"
)

ChIPseeker::plotPeakProf(
    peak = ranges,
    tagMatrix = tag_matrix,
    upstream = upstream, # rel(0.2),
    downstream = downstream, # rel(0.2),
    weightCol = "V5",
    ignore_strand = TRUE,
    conf = 0.95,
    by = "gene",
    type = "body",
    TxDb = txdb,
    facet = "row",
    nbin = 800,
    verbose = TRUE,
    resample = resample
)

dev.off()
base::message("Plot saved")

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()