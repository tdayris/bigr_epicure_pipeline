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
base::library(
    package = "TxDb.Mmusculus.UCSC.mm10.knownGene", character.only = TRUE
)
base::library(
    package = "TxDb.Hsapiens.UCSC.hg38.knownGene", character.only = TRUE
)
base::message("Libraries loaded")

tx_db <- TxDb.Hsapiens.UCSC.hg38.knownGene
if ("organism" %in% base::names(snakemake@params)) {
    if (base::as.character(x = snakemake@params[["organism"]]) == "mm10") {
        tx_db <- TxDb.Mmusculus.UCSC.mm10.knownGene
    }
}
de_counts <- base::readRDS(
    file = base::as.character(x = snakemake@input[["de_counts"]])
)

grouping_method <- "windows"
if ("grouping_method" %in% base::names(snakemake@params[["grouping_method"]])) {
    grouping_method <- bae::as.character(
        x = snakemake@params[["grouping_method"]]
    )
}

fdr <- NULL
if (grouping_method == "windows") {
    # Merging window method, aka: Quick and Dirty
    # Does not take clusters, annotations, nor interdependency
    # into account.
    base::message("Correcting multiple testing with binning")
    fdr <- csaw::mergeResults(
        ranges = de_counts,
        tab = rowData(de_counts),
        tol = base::as.integer(x = snakemake@params[["window_size"]])
    )

} else if (grouping_method == "genes") {
    # Merging method based on annotation
    base::message("Correcting multiple testing using genome annotation")

    broads <- genes(x = tx_db)
    broads <- resize(broads, width(broads) + 3000, fix = "end")

    fdr <- csaw::overlapResults(
        ranges = de_counts,
        regions = broads,
        tab = rowData(filtered.data)
    )

} else {
    base::stop("Unknown grouping method")
}

if ("rds" %in% base::names(snakemake@output)) {
    base::message("Saving RDS object")
    base::saveRDS(
        object = fdr,
        file = base::as.character(x = snakemake@output[["rds"]])
    )
}

if ("tsv" %in% base::names(snakemake@output)) {
    base::message("Saving results as TSV")
    utils::write.table(
        x = base::as.data.frame(fdr),
        file = base::as.character(x = snakemake@output[["tsv"]]),
        sep = "\t",
        quote = FALSE,
        col.names = TRUE,
        row.names = TRUE
    )
}


if ("csv" %in% base::names(snakemake@output)) {
    base::message("Saving results as CSV")
    utils::write.table(
        x = base::as.data.frame(fdr),
        file = base::as.character(x = snakemake@output[["csv"]]),
        sep = ",",
        quote = FALSE,
        col.names = TRUE,
        row.names = TRUE
    )
}
base::message("Process over")

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message");
base::sink();