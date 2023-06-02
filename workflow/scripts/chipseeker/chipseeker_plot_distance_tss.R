# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script filters counted reads with csaw

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")


# Load libraries
base::library(package = "ChIPseeker", character.only = TRUE)
base::message("Libraries loaded")

ranges <- readRDS(
    file = base::as.character(x = snakemake@input[["ranges"]])
)

# Build plot
png(
  filename = snakemake@output[["png"]],
  width = 1024,
  height = 768,
  units = "px",
  type = "cairo"
);

ChIPSeeker::plotDistToTSS(
    peakAnno,
    title = "Distribution of transcription factor-binding loci\nrelative to TSS"
)


dev.off()

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()