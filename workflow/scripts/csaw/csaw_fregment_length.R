# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script computes fragment length

base::library(package = "BiocParallel", character.only = TRUE)
base::library(package = "csaw", character.only = TRUE)
base::message("Library loaded")

# Loading input files
read_params <- base::readRDS(
    file = base::as.character(x = snakemake@input[["read_params"]])
)
read_params <- csaw::reform(x = read_params, dedup = TRUE)
base::message("Read parameters loaded, with dedup set to TRUE")

design <- base::readRDS(file = snakemake@input[["design"]])
base::message("Experimental design loaded with BamFiles")

max_delay <- 800

if (snakemake@threads > 1) {
    BiocParallel::register(
      BPPARAM = BiocParallel::MulticoreParam(workers = snakemake@threads)
    )
    options("mc.cores" = snakemake@threads)
    base::message("Process multi-threaded")
}

# Compute fragment size with cross-correlation
correlation <- csaw::correlateReads(
    bam.files = design$BamPath,
    max.dist = max_delay,
    cross = TRUE,
    param = read_params
)

frag_length <- csaw::maximizeCcf(correlation, ignore = 100)
base::message("Fragment length computed")

base::saveRDS(
    object = frag_length,
    file = base::as.character(x = snakemake@output[["fragment_length"]])
)
base::message("Fragment length saved")

grDevices::png(
    filename = base::as.character(x = snakemake@output[["png"]]),
    width = 1024,
    height = 788
)

graphics::plot(
    x = 0:max_delay,
    y = correlation,
    xlab = "Delay (bp)",
    ylab = "CCF",
    type = "l"
)
graphics::abline(v = frag_length, col = "red")
graphics::text(
    paste(frag_length, "bp"),
    x = frag_length,
    y = min(correlation),
    pos = "4",
    col = "red"
)

grDevices::dev.off()
base::message("Plot saved, process over")

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()
