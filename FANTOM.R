
library(CAGEfightR)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(ggvenn)

FANTOM <- read.table("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/HeLa_FANTOM5_enhancers.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

colnames(FANTOM)[1] <- "chr"
colnames(FANTOM)[2] <- "start"
colnames(FANTOM)[3] <- "end"

FANTOM_gr <- GRanges(
  seqnames = FANTOM$chr,
  ranges = IRanges(start = FANTOM$start + 1, end = FANTOM$end)
)

library(ggvenn)


library(rtracklayer)

# Import the chain file (adjust path to where you saved it)
chain <- import.chain("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/hg19ToHg38.over.chain")

# Liftover FANTOM from hg19 to hg38
FANTOM_gr_hg38 <- liftOver(FANTOM_gr, chain)
FANTOM_gr_hg38 <- unlist(FANTOM_gr_hg38)

export(FANTOM_gr_hg38, "HeLa_FANTOM5_enhancers_hg38.bed", format = "BED")