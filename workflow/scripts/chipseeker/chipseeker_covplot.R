# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script plots genome coverage using CHiPSeeker

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")
base::message("Logging defined")

# Load libraries
base::library(package = "ChIPseeker", character.only = TRUE)
base::message("Libraries loaded")

peak <- base::list(
  ChIPseeker::readPeakFile(
    peakfile = base::as.character(x = snakemake@input[["bed"]]),
    as = "GRanges"
  )
)
sname <- NULL
if ("sample" %in% base::names(x = snakemake@wildcards)) {
  sname <- c(base::as.character(x = snakemake@wildcards[["sample"]]))
} else if ("model_name" %in% base::names(x = snakemake@wildcards)) {
  sname <- c(base::as.character(x = snakemake@wildcards[["model_name"]]))
}
base::names(x = peak) <- sname
base::message("Ranges loaded")

# Build plot
png(
  filename = snakemake@output[["png"]],
  width = 1024,
  height = 768,
  units = "px",
  type = "cairo"
)

ChIPseeker::covplot(peak = peak, weightCol = "V5")


dev.off()
base::message("Plot saved")

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()