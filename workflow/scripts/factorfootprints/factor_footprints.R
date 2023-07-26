base::message("Loading libraries ... ")
suppressPackageStartupMessages(library("ATACseqQC"))
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg38.knownGene"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
suppressPackageStartupMessages(library("phastCons100way.UCSC.hg38"))
suppressPackageStartupMessages(library("MotifDb"))
suppressPackageStartupMessages(library("ChIPpeakAnno"))
suppressPackageStartupMessages(library("Rsamtools"))
base::message("Libraries loaded.")


base::message("Acquiering bam file...")
bamfile <- BamFile(
    file = base::as.character(x = snakemake@input[["bam"]]),
    index = base::as.character(x = snakemake@input[["bai"]])
)
name <- base::as.character(x = snakemake@params[["name"]])
base::message("BamFiles identified")

base::message("Loading bam file...")
bamdata <- readBamFile(
    bamFile = bamfile$path,
    bigFile = TRUE,
    asMates = TRUE,
    tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
)
base::message("Bam file loaded")

base::message("Shifting bam...")
shiftedBamfile <- base::as.character(x = snakemake@output[["shifted"]])
bamdata <- shiftGAlignmentsList(
    gal = bamdata,
    outbam = shiftedBamfile
)
base::message("Shift OK")

base::message("Acquiering motif...")
motif_name <- base::as.character(x = snakemake@params[["motif"]])
motif <- query(MotifDb, c(motif_name))
motif <- as.list(motif)
print(motif[[1]], digits = 2)
base::message("Motif retrieved.")

base::message("plot Footprints...")
genome <- Hsapiens
png(
    filename = snakemake@output[["png"]],
    width = 1024,
    height = 768,
    units = "px"
)
sigs <- factorFootprints(
    shiftedBamfile,
    pfm = motif[[1]],
    genome = genome,
    min.score = "90%",
    seqlev = c(1:23, "X", "Y"),
    upstream = 100,
    downstream = 100
)
dev.off()
base::message("Done.")

base::save.image(file = base::as.character(x = snakemake@output[["rda"]]))
base::message("Process over")