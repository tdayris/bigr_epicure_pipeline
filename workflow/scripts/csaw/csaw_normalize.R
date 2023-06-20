# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script normalize counted reads with csaw

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")


# Load libraries
base::library(package = "csaw", character.only = TRUE)
base::message("Libraries loaded")

norm_method <- "composition"
if ("norm_method" %in% base::names(snakemake@params)) {
    norm_method <- base::as.character(x = snakemake@params[["norm_method"]])
}


counts <- base::readRDS(
    file = base::as.character(x = snakemake@input[["counts"]])
)
base::message("Counts loaded")


if (norm_method == "composition") {
    # Highly enriched regions consume more sequencing resources and
    # thereby suppress the representation of other regions. Differences
    # in the magnitude of suppression between libraries can lead to
    # spurious DB calls. Scaling by library size fails to correct for this as
    # composition biases can still occur in libraries of the same size.
    base::message("Normalizing composition bias")

    binned <- base::readRDS(
        file = base::as.character(x = snakemake@input[["bins"]])
    )

    norm_factors <- csaw::normFactors(
        object = binned,
        se.out = counts
    )
} else if (norm_method == "efficiency") {
    # Efficiency biases in ChIP-seq data refer to fold changes in enrichment
    # that are introduced by variability in IP efficiencies between libraries.
    base::message("Normalizing efficiency bias")

    counts <- csaw::normFactors(
        object = binned,
        se.out = TRUE
    )
}


base::saveRDS(
    object = counts,
    file = base::as.character(x = snakemake@output[["rds"]])
)
base::message("Process over")

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()