# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script annotates differential peak coverage

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

# Acquiring annotation
extra <- 'dataset=c("hsapiens_gene_ensembl")'
if ("exta" %in% base::names(x = snakemake@params)) {
    extra <- base::as.character(x = snakemake@params[["extra"]])
}

command <- base::paste(
    "MEDIPS::MEDIPS.getAnnotation(",
    extra,
    ")"
)
base::message("Retrieving annotation with:")
base::print(command)

annotation <- base::eval(base::parse(text = command))
base::message("Annotation retrieved")


# Loading dataset that is to be annotated
edger_results <- base::readRDS(
    file = snakemake@input[["rds"]]
)
base::message("Peaks loaded")


# Annotate regions
annotation <- MEDIPS::MEDIPS.setAnnotation(
    regions = edger_results,
    annotation = annotation
)
base::message("Peaks annotated")

# Save results
base::saveRDS(
    object = annotation,
    file = snakemake@output[["rds"]]
)

utils::write.table(
    x = annotation,
    file = base::as.character(x = snakemake@output[["tsv"]]),
    sep = "\t"
)
base::message("Process over")