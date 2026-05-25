#clean up the data so the sequences matches the genome
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(GenomeInfoDb)

gr_plus <- import("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancer_supportdata/TLDR_plus.bedGraph")

#the unique 
Donotkeep_1 <- c(setdiff(seqlevels(gr_plus), seqlevels(seqinfo(BSgenome.Hsapiens.UCSC.hg38))))

gr_plus <- gr_plus[!seqnames(gr_plus) %in% Donotkeep_1]

# Export cleaned file
export(gr_plus, "/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancer_supportdata/TLDR_plus_clean2.bedGraph")


gr_minus <- import("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancer_supportdata/TLDR_minus.bedGraph")

#the unique 
Donotkeep_2 <- c(setdiff(seqlevels(gr_minus), seqlevels(seqinfo(BSgenome.Hsapiens.UCSC.hg38))))

gr_minus <- gr_minus[!seqnames(gr_minus) %in% Donotkeep_2]

# Export cleaned file
export(gr_minus, "/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancer_supportdata/TLDR_minus_clean2.bedGraph")
