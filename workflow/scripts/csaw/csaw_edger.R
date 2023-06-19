# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script search for significative differential peak calling
# with csaw/edgeR

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")


# Load libraries
base::library(package = "edgeR", character.only = TRUE)
base::library(package = "csaw", character.only = TRUE)
base::message("Libraries loaded")

counts <- base::readRDS(
    file = base::as.character(x = snakemake@input[["counts"]])
)
y <- csaw::asDGEList(object = counts)
base::message("Normalized-filtered counts loaded")

exp_design <- utils::read.table(
    file = base::as.character(x = snakemake@input[["design"]]),
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
)
design <- stats::model.matrix(
    object = stats::as.formula(snakemake@params[["formula"]]),
    data = exp_design
)
base::message("Model matrix built")

prior_n <- edgeR::getPriorN(
    y = y,
    design = design,
    prior.df = 20
)
base::message("Prior N estimated with 20 dof")

dispersion <- edgeR::estimateDisp(
    y = y,
    design = design,
    prior.n = prior_n
)
base::message("Dispersion estimated")

fit <- edgeR::glmQLFit(
    y = dispersion,
    design = design,
    robust = TRUE
)
base::message("GLM model fitted")

results <- edgeR::glmQLFTest(
    glmfit = fit,
    coef = utils::tail(x = base::colnames(x = design), 1)
)
base::message("Window counts tested")

rowData(counts) <- base::cbind(rowData(counts), results$table)
base::message("Results stored in RangedSummarizedExperiment object")

if ("csaw" %in% base::names(snakemake@output)) {
    base::saveRDS(
        object = counts,
        file = base::as.character(x = snakemake@output[["csaw"]])
    )
    base::message("csaw results saved")
}

if ("edger" %in% base::names(snakemake@ouptput)) {
    base::saveRDS(
        object = results,
        file = base::as.character(x = snakemake@output[["edger"]])
    )
    base::message("edgeR results saved")
}
base::message("Process over")

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message");
base::sink();