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
base::library(package = "DNAcopy", character.only = TRUE)
base::message("Libraries loaded")


mset1 <- base::readRDS(file = base::as.character(x = "mset1"))
extra <- "MSet1 = mset1, diff.method = 'edgeR', p.adj = 'bonferroni'"
    base::message("MSet1 loaded")

mset2 <- NULL
if ("mset2" %in% base::names(snakemake@input)) {
    mset2 <- base::readRDS(file = base::as.character(x = "mset2"))
    extra <- base::paste(extra, "MSet2 = mset2", sep = ", ")
    base::message("MSet2 loaded")
}

cset <- NULL
if ("cset" %in% base::names(snakemake@input)) {
    cset <- base::readRDS(file = base::as.character(x = "cset"))
    extra <- base::paste(extra, "CSet = cset", sep = ", ")
    base::message("CSet loaded")
}

iset1 <- NULL
if ("iset1" %in% base::names(snakemake@input)) {
    iset1 <- iset1base::readRDS(file = base::as.character(x = "iset1"))
    extra <- base::paste(extra, "ISet1 = iset1", sep = ", ")
    base::message("ISet1 loaded")
}

iset2 <- NULL
if ("iset2" %in% base::names(snakemake@input)) {
    iset2 <- iset1base::readRDS(file = base::as.character(x = "iset2"))
    extra <- base::paste(extra, "ISet2 = iset2", sep = ", ")
    base::message("ISet2 loaded")
}

if (! (base::is.null(iset1) && base::is.null(iset2))) {
    extra <- base::paste(extra, "CNV = TRUE", sep = ", ")
    base::message(
        "copy number variation will be tested by ",
        "applying the package DNAcopy to the ",
        "window-wise log-ratios calculated based on ",
        "the the means per group"
    )
}

if (! (base::is.null(mset2) && base::is.null(cset))) {
    extra <- base::paste(extra, "MeDIP = TRUE", sep = ", ")
    base::message(
        "CpG density dependent relative methylation scores ",
        "(rms) will be calculated for the MEDIPS SETs"
    )
}

medip_meth_command <- base::paste0(
    "MEDIPS::MEDIPS.meth(",
    extra,
    ")"
)
base::message("Command line used to determine coverage and analysis methods:")
base::message(medip_meth_command)

coverage_profiles <- base::eval(base::parse(text = medip_meth_command))

base::saveRDS(
    object = coverage_profiles,
    file = base::as.character(x = snakemake@output[["rds"]])
)

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
sink(type = "message")
sink()