base::message("Loading libraries ... ")
suppressPackageStartupMessages(library("ATACseqQC"))
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg38.knownGene"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
suppressPackageStartupMessages(library("phastCons100way.UCSC.hg38"))
suppressPackageStartupMessages(library("MotifDb"))
suppressPackageStartupMessages(library("ChIPpeakAnno"))
suppressPackageStartupMessages(library("Rsamtools"))
base::message("Libraries loaded.")

# base::message("Setting sequence level style...")
# seqlevelsStyle(TxDb.Hsapiens.UCSC.hg38.knownGene) <- "Ensembl"
# seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "Ensembl"
# base::message("Database chromosome renamed.")

base::message("Acquiering bam file...")
bamfile <- BamFile(
    file = base::as.character(x = snakemake@input[["bam"]])
)
name <- base::as.character(x = snakemake@params[["name"]])
base::print(bamfile)
base::print(name)
base::message("BamFiles identified")

base::message("Reading bam tags...")
possibleTag <- list(
    "integer" = c(
        "AM", "AS", "CM", "CP", "FI", "H0", "H1",
        "H2", "HI", "IH", "MQ", "NH", "NM", "OP",
        "PQ", "SM", "TC", "UQ"
    ),
    "character" = c(
        "BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
        "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
        "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
        "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS", "U2"
    )
)
bamTop100 <- scanBam(
    BamFile(bamfile$path, yieldSize = 100),
    param = ScanBamParam(tag=unlist(possibleTag))
)[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100) > 0]
base::print(tags)
base::message("Tags Acquired")

base::message("Retrieving sequence level informations...")
seqlev <- as.vector(
    sapply(c(1:22, "X", "Y"), function(chrom) paste0("chr", chrom))
)
seqinformation <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
which <- as(seqinformation[seqlev], "GRanges")
base::print(which)
base::message("Sequences retrived")

base::message("Loading bam file...")
bamdata <- readBamFile(
    bamFile = bamfile$path,
    bigFile = TRUE,
    asMates = TRUE,
    tags = tags,
    which = which,
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
    seqlev = seqlev,
    upstream = 100,
    downstream = 100
)
dev.off()
base::message("Done.")

base::save.image(file = base::as.character(x = snakemake@output[["rda"]]))
base::message("Process over")
