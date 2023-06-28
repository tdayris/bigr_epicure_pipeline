# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script plots genomic annotation using CHiPSeeker
# on a list of ranges

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")
base::message("Logging defined")

# Load libraries
base::library(package = "ChIPseeker", character.only = TRUE)
base::message("Libraries loaded")

ranges_list <- base::lapply(
    snakemake@input[["ranges"]],
    function(range_path) {
        base::readRDS(
            file = base::as.character(x = range_path)
        )
    }
)
base::message("Ranges loaded")

# Build plot
grDevices::png(
  filename = snakemake@output[["png"]],
  width = 1024,
  height = 384,
  units = "px",
  type = "cairo"
)

ChIPseeker::plotAnnoBar(ranges_list)


dev.off()
base::message("Plot saved")

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()