# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script loads bam file into a MEDIPS dataset

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")


# Load libraries
base::library(package = "MEDIPS", character.only = TRUE)
base::library(package = "BSgenome.Hsapiens.UCSC.hg38", character.only = TRUE)
base::library(package = "BSgenome.Mmusculus.UCSC.mm10", character.only = TRUE)
base::message("Libraries loaded")

medips_sets <- c()

for (s in seq_along(along.with = snakemake@input[["bam"]])) {
    # Load optional parameters
    bam_path <- base::as.character(x = snakemake@input[["bam"]][s])

    medips_extra <- 'file = bam_path'

    if ("extra" %in% base::names(snakemake@params)) {
        medips_extra <- base::paste(
            medips_extra,
            base::as.character(x = snakemake@params[["extra"]][s]),
            sep = ", "
        )
    }

    medips_createset_command <- base::paste0(
        "MEDIPS::createSet(",
        medips_extra,
        ")"
    )
    base::message("Command line used to load ", snakemake@input[["bam"]][s])
    base::message(medips_createset_command)

    base::append(
        x = medips_sets,
        values = base::eval(base::parse(text = medips_createset_command))
    )
    base::message("MEDIPS-set created")
}

base::saveRDS(
    x = medips_sets,
    file = base::as.character(x = snakemake@output[["rds"]])
)
base::message("Process over")


# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
sink(type = "message")
sink()