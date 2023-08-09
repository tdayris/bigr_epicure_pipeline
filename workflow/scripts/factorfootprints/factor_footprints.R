base::message("Loading libraries ... ")
suppressPackageStartupMessages(library("ATACseqQC"))
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg38.knownGene"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
suppressPackageStartupMessages(library("phastCons100way.UCSC.hg38"))
suppressPackageStartupMessages(library("MotifDb"))
suppressPackageStartupMessages(library("ChIPpeakAnno"))
suppressPackageStartupMessages(library("Rsamtools"))
base::message("Libraries loaded.")

base::message("Setting sequence level style...")
seqlevelsStyle(TxDb.Hsapiens.UCSC.hg38.knownGene) <- "Ensembl"
seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "Ensembl"
base::message("Database chromosome renamed.")

base::message("Acquiering bam file...")
bamfile <- BamFile(
    file = base::as.character(x = snakemake@input[["bam"]]),
    index = base::as.character(x = snakemake@input[["bai"]])
)
name <- base::as.character(x = snakemake@params[["name"]])
print(bamfile)
print(name)
base::message("BamFiles identified")

base::message("Loading bam file...")
bamdata <- readBamFile(
    bamFile = bamfile$path,
#    index = bamfile$path,
    bigFile = TRUE,
    asMates = TRUE,
    tags = c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
)
base::message("Bam file loaded")

base::message("Shifting bam...")
shiftedBamfile <- base::as.character(x = snakemake@output[["shifted"]])
shiftedBamdir <- base::dirname(shiftedBamfile)
print(shiftedBamdir)
base::dir.create(
    path = shiftedBamdir,
    recursive = TRUE
)
bamdata <- shiftGAlignmentsList(
    gal = bamdata,
    outbam = shiftedBamfile
)
print(bamdata)
base::message("Shift OK")

base::message("Acquiering motif...")
motif_name <- base::as.character(x = snakemake@params[["motif"]])
motif <- query(MotifDb, c(motif_name))
motif <- as.list(motif)
print(motif[[1]], digits = 2)
base::message("Motif retrieved.")

base::message("plot Footprints...")
genome <- Hsapiens
print(genome)

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
    seqlev = c(1:22, "X", "Y"),
    upstream = 100,
    downstream = 100
)
dev.off()
base::message("Done.")

base::save.image(file = base::as.character(x = snakemake@output[["rda"]]))
base::message("Process over")
