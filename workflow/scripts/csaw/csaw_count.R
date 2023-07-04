# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script counts reads over sliding windows with csaw

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")

# Load libraries
base::library(package = "BiocParallel", character.only = TRUE)
base::library(package = "csaw", character.only = TRUE)
base::message("Libraries loaded")

if (snakemake@threads > 1) {
    BiocParallel::register(
        BiocParallel::MulticoreParam(snakemake@threads)
    )
    options("mc.cores" = snakemake@threads)
    base::message("Process multi-threaded")
}

# Loading list of input files
design <- readRDS(
    file = base::as.character(x = snakemake@input[["design"]])
)

frag_length <- readRDS(
    file = base::as.character(x = snakemake@input[["fragment_length"]])
)

# Loading filter parameters
read_params <- base::readRDS(
    file = base::as.character(x = snakemake@input[["read_params"]])
)
base::message("Input data loaded")

extra <- "bam.files = design$BamPath, param = read_params, ext=frag.len"
if ("extra" %in% base::names(x = snakemake@params)) {
    extra <- base::paste(
        extra,
        base::as.character(x = snakemake@params[["extra"]]),
        sep = ", "
    )
}


command <- base::paste(
    "csaw::windowCounts(",
    extra,
    ")"
)
base::message("Count command line:")
base::print(command)

# Counting reads over sliding window
counts <- base::eval(base::parse(text = command))
base::message(
    "Number of extended reads ",
    "overlapping a sliding window ",
    "at spaced positions across the ",
    "genome, acquired."
)

# Saving results
base::saveRDS(
    object = counts,
    file = snakemake@output[["rds"]]
)
base::print(counts$totals)
base::message("RDS saved, process over")


# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()