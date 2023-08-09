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
base::library(package = "csaw", character.only = TRUE)
base::message("Libraries loaded")


counts <- base::readRDS(
    file = base::as.character(x = snakemake@input[["counts"]])
)

log_threshold <- 1.1
if ("log_threshold" %in% base::names(snakemake@params)) {
    log_threshold <- as.numeric(snakemake@params[["log_threshold"]])
}


filter_method <- "average_log_cpm"
if ("filter_method" %in% base::names(x = snakemake@params)) {
    filter_method <- base::as.character(
        x = snakemake@params[["filter_method"]]
    )
}

keep <- base::data.frame(
    "Kept" = 0,
    "Filtered" = 0,
    "PercentFiltered" = 0,
    "PercentKept" = 0
)
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
    keep <- abundances > aveLogCPM(min_cpm, lib.size = mean(counts$totals))
    counts <- counts[keep, ]
} else if (filter_method == "proportion") {
    # Filter based on estimated proportion of signal of interest / noise
    # Not that good if site proportion varies across samples

    base::message("Filtering based on signal / noise proportion")

    proportion <- 0.999
    if ("proportion" %in% base::names(x = snakemake@params)) {
        proportion <- base::as.numeric(x = snakemake@params[["proportion"]])
    }

    keep <- csaw::filterWindowsProportion(counts)$filter > proportion
    counts <- counts[keep, ]
} else if (filter_method == "global_enrichment") {
    # Filtering by global enrichment. The degree of background enrichment
    # is estimated by counting reads into large bins across the genome.
    # Not that good if there are differences in the background coverage
    # across the genome.

    base::message("Filtering based on global enrichment")

    binned <- base::readRDS(
        file = base::as.character(x = snakemake@input[["binned"]])
    )

    filter_stat <- csaw::filterWindowsGlobal(
        data = counts,
        background = binned
    )

    keep <- filter_stat$filter > log2(log_threshold)
    counts <- counts[keep, ]
} else if (filter_method == "local_enrichment") {
    # Mimicking single-sample peak callers.
    # Not good if any peaks are in the neighborhood ! Please make
    # sure there are not peak clusters of any kind.

    base::message("Filtering based on local enrichment")

    neighbor <- base::readRDS(
        file = base::as.character(x = snakemake@input[["neighbor"]])
    )

    filter_stat <- csaw::filterWindowsLocal(
        data = counts,
        background = neighbor
    )

    keep <- filter_stat$filter > log2(log_threshold)
    counts <- counts[keep, ]
} else if (filter_method == "local_maxima") {
    # Use highest average abundance within 1kbp on either side.
    # Very aggressive, use only in datasets with sharp binding.

    base::message("Filtering based on local peak maxima")

    maxed <- base::readRDS(
        file = base::as.character(x = snakemake@input[["maxed"]])
    )

    filter_stat <- csaw::filterWindowsLocal(
        data = counts,
        background = maxed
    )

    keep <- filter_stat$filter > log2(log_threshold)
    counts <- counts[keep, ]
} else if (filter_method == "input_controls") {
    # Use negative controls for ChIP-seq to account for
    # both sequencing and mapping biases in ChIP-seq data

    base::message("Filtering based on negative controls (aka Input)")

    input_counts <- base::readRDS(
        file = base::as.character(x = snakemake@input[["input_counts"]])
    )
    base::message("Input counts loaded")
    counts_binned <- base::readRDS(
        file = base::as.character(x = snakemake@input[["binned"]])
    )
    base::message("Tested (binned) counts loaded")
    input_binned <- base::readRDS(
        file = base::as.character(x = snakemake@input[["input_binned"]])
    )
    base::message("Input (binned) counts loaded")

    scale_info <- csaw::scaleControlFilter(
        data.bin = counts_binned,
        back.bin = input_binned
    )
    base::message("Scale control filter acquired")

    filter_stat <- csaw::filterWindowsControl(
        data = counts,
        background = input_counts,
        prior.count = 2,
        scale.info = scale_info
    )
    base::message("Negative control signal filtered")

    keep <- filter_stat$filter > log2(log_threshold)
    counts <- counts[keep, ]
} else {
    base::stop("Unknown filter method")
}


k <- base::summary(object = keep)
keep <- base::as.data.frame()
keep$Kept <- k[["TRUE"]]
keep$Filtered <- k[["FALSE"]]
keep$PercentFiltered <- keep$Filtered / (keep$Kept + keep$Filtered)
keep$PercentKept <- keep$Kept / (keep$Kept + keep$Filtered)


# Saving results
utils::write.table(
    x = keep,
    file = base::as.character(x = snakemake@output[["qc"]]),
    sep = "\t"
)
base::message("QC table saved")

if ("png" %in% base::names(x = snakemake@output)) {
    # Saving QC plot
    grDevices::png(
        filename = base::as.character(x = snakemake@output[["png"]]),
        width = 1024,
        height = 768
    )

    graphics::hist(
        filter_stat$filter,
        main = "",
        breaks = 50,
        xlab = "Background abundance (log2-CPM)"
    )
    graphics::abline(v = log2(log_threshold), col = "red")

    grDevices::dev.off()
    base::message("QC plot saved")
}


base::saveRDS(
    object = counts,
    file = base::as.character(x = snakemake@output[["rds"]])
)
base::message("RDS saved, process over")

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()