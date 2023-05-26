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
base::library(package = "edgeR", character.only = TRUE)
base::library(package = "csaw", character.only = TRUE)
base::message("Libraries loaded")


counts <- base::readRDS(file = base::as.character(x = snakemake@input[["rds"]]))


filter_method <- "average_log_cpm"
if ("filter_method" %in% base::names(x = snakemake@params)) {
    filter_method <- base::as.character(x = snakemake@params[["filter_method"]])
}
base::message("Input data loaded")

# Filtering counts
if (filter_method == "average_log_cpm") {
    # Filter based on average CPM counts.
    # not that good in case of large sequencing libraries
    # not that good in case of large variability in sequencing library sizes

    base::message("Filtering based on count size (aveLogCPM)")
    min_cpm <- 5
    if ("min_cpm" %in% base::names(snakemake@params)) {
        min_cpm <- base::as.numeric(x = snakemake@params[["min_cpm"]])
    }
    abundance <- edgeR::aveLogCPM(SummarizedExperiment::asDGEList(counts))
    keep <- abundances > aveLogCPM(min_cpm, lib.size=mean(counts$totals))
    counts <- counts[keep,]
} else if (filter_method == "proportion") {
    # Filter based on estimated proportion of signal of interest / noise
    # Not that good if site proportion varies across samples
    proportion <- 0.999
    if ("proportion" %in% base::names(x = snakemake@params)) {
        proportion <- base::as.numeric(x = snakemake@params[["proportion"]])
    }

    keep <- csaw::filterWindowsProportion(counts)$filter > proportion
    counts <- counts[keep,]
}

# Saving results
base::saveRDS(
    object = counts, 
    file = base::as.character(x = snakemake@output[["rds"]])
)

