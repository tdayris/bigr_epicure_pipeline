# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script builds ReadParam object for csaw

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")

# Load libraries
base::library(package = "csaw", character.only = TRUE)
base::library(package = "rtracklayer", character.only = TRUE)
base::message("Libraries loaded")

standard_chrom <- c(1:23, "X", "Y")
if ("organism" %in% base::names(x = snakemake@params)) {
    if (base::as.character(x = snakemake@params[["organism"]]) == "mm10") {
        standard_chrom <- c(1:19, "X", "Y")
    }
}

blacklist <- base::as.character(x = snakemake@input[["blacklist"]])
blacklist <- rtracklayer::import(blacklist)

# Build command line
extra <- "minq=20, restrict=standard_chrom, discard=blacklist"
if ("extra" %in% base::names(x = snakemake@params)) {
    extra <- base::as.character(x = snakemake@params[["extra"]])
}

command <- base::paste0(
    "csaw::readParam(",
    extra,
    ")"
)
base::message("Command line:")
base::print(command)


# Run command line
read_params <- base::eval(base::parse(text = command))

# Save results
base::saveRDS(
    object = read_params,
    file = snakemake@output[["rds"]]
)
base::message("RDS saved, process over")

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()