# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script annotates peaks with CHiPSeeker

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")


# Load libraries
base::library(package = "ChIPseeker", character.only = TRUE)

base::library(package = "org.Hs.eg.db", character.only = TRUE)
base::library(package = "org.Mm.eg.db", character.only = TRUE)

base::library(package = "EnsDb.Hsapiens.v86", character.only = TRUE)
base::library(package = "EnsDb.Mmusculus.v79", character.only = TRUE)
base::message("Libraries loaded")

organism <- "hg38"
if ("organism" %in% base::names(x = snakemake@params)) {
    organism <- base::as.character(x = snakemake@params[["organism"]])
}

organism <- base::tolower(x = organism)
if (organism == "hg38") {
    anno_db <- "org.Hs.eg.db"
    edb <- EnsDb.Hsapiens.v86
} else if (organism == "mm10") {
    anno_db <- "org.Mm.eg.db"
    edb <- EnsDb.Mmusculus.v79
} else {
    base::stop("Unknown organism annotation")
}
seqlevelsStyle(edb) <- "Ensembl"

ranges <- NULL
if ("ranges" %in% base::names(x = snakemake@input)) {
    ranges <- base::readRDS(file = base::as.character(x = snakemake@input[["ranges"]]))
} else {
    ranges <- ChIPseeker::readPeakFile(
        peakfile = base::as.character(x = snakemake@input[["bed"]]),
        as = "GRanges"
    )
}
base::message("Parameters and data loaded")


# Building command line
extra <- "peak = ranges, annoDb = anno_db, TxDb = edb"
if ("extra" %in% base::names(x = snakemake@params)) {
    extra <- base::paste(
        extra,
        base::as.character(x = snakemake@params[["extra"]]),
        sep = ", "
    )
}

command <- base::paste0(
    "ChIPseeker::annotatePeak(",
    extra,
    ")"
)
base::message("Command line used:")
base::message(command)

# Annotating
annotation <- base::eval(base::parse(text = command))
base::message("Data annotated")


# Saving results
base::saveRDS(
    object = annotation,
    file = base::as.character(x = snakemake@output[["rds"]])
)
base::message("RDS saved")

utils::write.table(
    base::as.data.frame(annotation),
    sep = "\t",
    col.names = TRUE,
    row.names = TRUE,
    file = base::as.character(x = snakemake@output["tsv"])
)
base::message("Provess over")


# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()