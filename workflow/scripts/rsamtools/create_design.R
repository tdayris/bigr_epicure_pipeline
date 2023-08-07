# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script load several bamfiles and
# builds a design file for future use

base::library(package = "Rsamtools", character.only = TRUE)
base::message("Library loaded")

# Loading input data
rds_paths <- base::sapply(
    snakemake@input[["bam_files"]],
    function(rds_path) base::readRDS(file = base::as.character(x = rds_path))
)
base::message("BamFiles loaded")

design <- utils::read.table(
    file = base::as.character(x = snakemake@input[["design"]]),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)
base::rownames(design) <- design$Sample_id
base::message("Design loaded")

# Building quality control table
design$BamPath <- rds_paths
base::message("Design enhanced")

diagnostics <- base::list()
for (i in base::seq_along(design$BamPath)) {
    bam <- design$BamPath[[i]]
    total_reads <- Rsamtools::countBam(file = bam)$records

    params <- Rsamtools::ScanBamParam(
        flag = Rsamtools::scanBamFlag(isUnmapped = FALSE)
    )
    tota_mapped <- Rsamtools::countBam(file = bam, param = params)

    params <- Rsamtools::ScanBamParam(
        flag = Rsamtools::scanBamFlag(isUnmapped = FALSE, isDuplicate = TRUE)
    )
    marked <- Rsamtools::countBam(bam, param = params)$records

    diagnostics[[i]] <- c(
        Total = total_reads, 
        Mapped = tota_mapped, 
        Marked = marked
    )
}

base::print(diagnostics)
diag_stats <- base::data.frame(base::do.call(rbind, diagnostics))
base::rownames(diag_stats) <- design$Sample_id


diag_stats$Prop.mapped <- (as.numeric(x = diag_stats[["Mapped.records"]]) / as.numeric(x = diag_stats[["Total"]])) * 100
diag_stats$Prop.marked <- (as.numeric(x = diag_stats[["Marked"]]) / as.numeric(x = diag_stats[["Mapped.records"]])) * 100
base::message("Quality controls performed")

# diag_stats <- data.frame(t(sapply(diag_stats, c)))
diag_stats <- data.frame(lapply(diag_stats, as.character), stringsAsFactors=FALSE)
qc_path <- base::as.character(x = snakemake@output[['qc']])
print(diag_stats)
print(qc_path)
print(str(diag_stats))
# Saving results
utils::write.table(
    x = diag_stats,
    file = qc_path,
    sep = "\t",
    row.names = TRUE,
    col.names = TRUE
)
base::message("QC saved")

base::saveRDS(
    object = design,
    file = base::as.character(x = snakemake@output[["design"]])
)
base::message("Design saved")

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()
