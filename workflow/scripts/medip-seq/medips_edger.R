# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script computes differential peak coverage over MEDIPS dataset

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")


# Load libraries
base::library(package = "MEDIPS", character.only = TRUE)
base::library(package = "BSgenome.Hsapiens.UCSC.hg38", character.only = TRUE)
base::library(package = "BSgenome.Mmusculus.UCSC.mm10", character.only = TRUE)

results <- base::readRDS(
    file = base::as.character(x = snakemake@input[["meth"]])
)
base::message("Dataset loaded")


# Build command line
extra <- "results = results"
if ("extra" %in% base::names(x = snakemake@params)) {
    extra <- base::paste(
        extra,
        base::as.character(x = snakemake@params[["extra"]]),
        sep = ", "
    )
}

command <- base::paste0(
    "MEDIPS.selectSig(",
    extra,
    ")"
)
base::message("Command line:")
base::paste(command)


# Launch analysis
edger_results <- base::eval(base::parse(text = command))
base::message("Differential peak coverage perfomed")

# Optionally merge close peaks
if ("merge_distance" %in% base::names(x = snakemake@params)) {
    edger_results <- MEDIPS.mergeFrames(
        frames = edger_results,
        distance = base::as.numeric(x = snakemake@params[["merge_distance"]])
    )
    base::message("Close peaks merged")
}

# Save results
base::saveRDS(
    object = edger_results,
    file = base::as.character(x = snakemake@output[["rds"]])
)

utils::write.table(
    x = edger_results,
    file = base::as.character(x = snakemake@output[["tv"]]),
    sep = "\t",
    quote = FALSE
)
base::message("Process over")