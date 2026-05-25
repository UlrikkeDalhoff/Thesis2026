library(eulerr)

makeVenn <- function(mylist)
{
  all <- unique(unlist(mylist))
  all <- all[!(all=="")]
  dmat <- matrix(0,ncol=length(mylist),nrow=length(all))
  row.names(dmat) <- all
  colnames(dmat) <- names(mylist)
  for(i in 1:length(mylist))
  {
    mgenes <- mylist[[i]]
    mgenes <- mgenes[!(mgenes=="")]
    dmat[mgenes,i] <- 1
  }
  genes.venn <- euler(dmat, shape = "ellipse")  
  plot(genes.venn, quantities = TRUE)
}




library(CAGEfightR)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(ggvenn)


#We create a CTSS object and get a count at each 5' position
ctss_exp1 <- quantifyCTSSs(
  plusStrand = c(sample1 = "/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancers_from_ourdata/plus.bedGraph"),
  minusStrand = c(sample1 = "/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancers_from_ourdata/minus.bedGraph"),
  genome = seqinfo(BSgenome.Hsapiens.UCSC.hg38)
)

#We create a CTSS object and get a count at each 5' position
ctss_exp2 <- quantifyCTSSs(
  plusStrand = c(sample1 = "/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancer_supportdata/TLDR_plus_clean2.bedGraph"),
  minusStrand = c(sample1 = "/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancer_supportdata/TLDR_minus_clean2.bedGraph"),
  genome = seqinfo(BSgenome.Hsapiens.UCSC.hg38)
)

ctss_TLDR <- quantifyCTSSs(
  plusStrand = c(sample1 = "/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancer_supportdata/TLDR_plus_clean2.bedGraph"),
  minusStrand = c(sample1 = "/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancer_supportdata/TLDR_minus_clean2.bedGraph"),
  genome = seqinfo(BSgenome.Hsapiens.UCSC.hg38)
)

ctss_SLIC <- quantifyCTSSs(
  plusStrand = c(sample1 = "/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancer_supportdata/SLIC_plus.bedGraph"),
  minusStrand = c(sample1 = "/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancer_supportdata/SLIC_minus.bedGraph"),
  genome = seqinfo(BSgenome.Hsapiens.UCSC.hg38)
)


#And then we normalize
ctss_exp1 <- calcTPM(ctss_exp1)

#Combine samples 
ctss_exp1 <- calcPooled(ctss_exp1)

#And then we normalize
ctss_exp2 <- calcTPM(ctss_exp2)

#Combine samples 
ctss_exp2 <- calcPooled(ctss_exp2)

#Normalize 
ctss_TLDR <- calcTPM(ctss_TLDR)

#Combine
ctss_TLDR <- calcPooled(ctss_TLDR)

#Normalize 
ctss_SLIC <- calcTPM(ctss_SLIC)

#Combine
ctss_SLIC <- calcPooled(ctss_SLIC)




#Now we can detect bidirectional clusteres (enhancers)
bc075_100_exp1 <- clusterBidirectionally(
  ctss_exp1,
  balanceThreshold = 0.75,
  window = 101
)

bc075_200_exp1 <- clusterBidirectionally(
  ctss_exp1,
  balanceThreshold = 0.75,
  window = 201
)

bc075_400_exp1 <- clusterBidirectionally(
  ctss_exp1,
  balanceThreshold = 0.75,
  window = 401
)

bc075_600_exp1 <- clusterBidirectionally(
  ctss_exp1,
  balanceThreshold = 0.75,
  window = 601
)

bc075_800_exp1 <- clusterBidirectionally(
  ctss_exp1,
  balanceThreshold = 0.75,
  window = 801
)


#Now we can detect bidirectional clusteres (enhancers)
bc075_100_exp2 <- clusterBidirectionally(
  ctss_exp2,
  balanceThreshold = 0.75,
  window = 101
)

bc075_200_exp2 <- clusterBidirectionally(
  ctss_exp2,
  balanceThreshold = 0.75,
  window = 201
)

bc075_400_exp2 <- clusterBidirectionally(
  ctss_exp2,
  balanceThreshold = 0.75,
  window = 401
)

bc075_600_exp2 <- clusterBidirectionally(
  ctss_exp2,
  balanceThreshold = 0.75,
  window = 601
)

bc075_800_exp2 <- clusterBidirectionally(
  ctss_exp2,
  balanceThreshold = 0.75,
  window = 801
)

#Now we can detect bidirectional clusteres (enhancers)
bc075_100_TLDR <- clusterBidirectionally(
  ctss_TLDR,
  balanceThreshold = 0.75,
  window = 101
)

bc075_200_TLDR <- clusterBidirectionally(
  ctss_TLDR,
  balanceThreshold = 0.75,
  window = 201
)

bc075_400_TLDR <- clusterBidirectionally(
  ctss_TLDR,
  balanceThreshold = 0.75,
  window = 401
)

bc075_600_TLDR <- clusterBidirectionally(
  ctss_TLDR,
  balanceThreshold = 0.75,
  window = 601
)

bc075_800_TLDR <- clusterBidirectionally(
  ctss_TLDR,
  balanceThreshold = 0.75,
  window = 801
)


#Now we can detect bidirectional clusteres (enhancers)
bc075_100_SLIC <- clusterBidirectionally(
  ctss_SLIC,
  balanceThreshold = 0.75,
  window = 101
)

bc075_200_SLIC <- clusterBidirectionally(
  ctss_SLIC,
  balanceThreshold = 0.75,
  window = 201
)

bc075_400_SLIC <- clusterBidirectionally(
  ctss_SLIC,
  balanceThreshold = 0.75,
  window = 401
)

bc075_600_SLIC <- clusterBidirectionally(
  ctss_SLIC,
  balanceThreshold = 0.75,
  window = 601
)

bc075_800_SLIC <- clusterBidirectionally(
  ctss_SLIC,
  balanceThreshold = 0.75,
  window = 801
)


library(GenomicRanges)
library(eulerr)


# Combine all regions
all_regions <- reduce(c(
  bc075_100_exp1,
  bc075_100_exp2,
  bc075_100_TLDR,
  bc075_100_SLIC
))

# Function to check overlap membership
get_membership <- function(query, subject) {
  hits <- findOverlaps(subject, query, ignore.strand = TRUE)
  unique(queryHits(hits))
}

# Membership indices in common universe 
exp1_idx <- get_membership(bc075_100_exp1, all_regions)
exp2_idx <- get_membership(bc075_100_exp2, all_regions)
TLDR_idx <- get_membership(bc075_100_TLDR, all_regions)
SLIC_idx <- get_membership(bc075_100_SLIC, all_regions)

# Convert to character
makeVenn(list(
  "Exp1 BT075 W100" = as.character(exp1_idx),
  "Exp2 BT075 W100" = as.character(exp2_idx),
  "TLDR BT075 W100" = as.character(TLDR_idx),
  "SLIC BT075 W100" = as.character(SLIC_idx)
))

#for BT075 window 200
# Combine all regions
all_regions_BT075_200 <- reduce(c(
  bc075_200_exp1,
  bc075_200_exp2,
  bc075_200_TLDR,
  bc075_200_SLIC
))


# Membership indices in common universe 
exp1_idx_BT075_200 <- get_membership(bc075_200_exp1, all_regions_BT075_200)
exp2_idx_BT075_200 <- get_membership(bc075_200_exp2, all_regions_BT075_200)
TLDR_idx_BT075_200 <- get_membership(bc075_200_TLDR, all_regions_BT075_200)
SLIC_idx_BT075_200 <- get_membership(bc075_200_SLIC, all_regions_BT075_200)

# Convert to character
makeVenn(list(
  "Exp1 BT075 W200" = as.character(exp1_idx_BT075_200),
  "Exp2 BT075 W200" = as.character(exp2_idx_BT075_200),
  "TLDR BT075 W200" = as.character(TLDR_idx_BT075_200),
  "SLIC BT075 W200" = as.character(SLIC_idx_BT075_200)
))

#for BT075 window 400
# Combine all regions
all_regions_BT075_400 <- reduce(c(
  bc075_400_exp1,
  bc075_400_exp2,
  bc075_400_TLDR,
  bc075_400_SLIC
))


# Membership indices in common universe 
exp1_idx_BT075_400 <- get_membership(bc075_400_exp1, all_regions_BT075_400)
exp2_idx_BT075_400 <- get_membership(bc075_400_exp2, all_regions_BT075_400)
TLDR_idx_BT075_400 <- get_membership(bc075_400_TLDR, all_regions_BT075_400)
SLIC_idx_BT075_400 <- get_membership(bc075_400_SLIC, all_regions_BT075_400)

# Convert to character
makeVenn(list(
  "Exp1 BT075 W400" = as.character(exp1_idx_BT075_400),
  "Exp2 BT075 W400" = as.character(exp2_idx_BT075_400),
  "TLDR BT075 W400" = as.character(TLDR_idx_BT075_400),
  "SLIC BT075 W400" = as.character(SLIC_idx_BT075_400)
))


makeVenn2 <- function(mylist)
{
  mylist <- lapply(mylist, function(x) x[x != ""])
  ggvenn(mylist, show_percentage = FALSE)
}

makeVenn2(list(
  "Exp1 BT075 W400" = as.character(exp1_idx_BT075_400),
  "Exp2 BT075 W400" = as.character(exp2_idx_BT075_400),
  "TLDR BT075 W400" = as.character(TLDR_idx_BT075_400),
  "SLIC BT075 W400" = as.character(SLIC_idx_BT075_400)
))


#for BT075 window 600
# Combine all regions
all_regions_BT075_600 <- reduce(c(
  bc075_600_exp1,
  bc075_600_exp2,
  bc075_600_TLDR,
  bc075_600_SLIC
))


# Membership indices in common universe 
exp1_idx_BT075_600 <- get_membership(bc075_600_exp1, all_regions_BT075_600)
exp2_idx_BT075_600 <- get_membership(bc075_600_exp2, all_regions_BT075_600)
TLDR_idx_BT075_600 <- get_membership(bc075_600_TLDR, all_regions_BT075_600)
SLIC_idx_BT075_600 <- get_membership(bc075_600_SLIC, all_regions_BT075_600)

# Convert to character
makeVenn(list(
  "Exp1 BT075 W600" = as.character(exp1_idx_BT075_600),
  "Exp2 BT075 W600" = as.character(exp2_idx_BT075_600),
  "TLDR BT075 W600" = as.character(TLDR_idx_BT075_600),
  "SLIC BT075 W600" = as.character(SLIC_idx_BT075_600)
))



#for BT075 window 800
# Combine all regions
all_regions_BT075_800 <- reduce(c(
  bc075_800_exp1,
  bc075_800_exp2,
  bc075_800_TLDR,
  bc075_800_SLIC
))


# Membership indices in common universe 
exp1_idx_BT075_800 <- get_membership(bc075_800_exp1, all_regions_BT075_800)
exp2_idx_BT075_800 <- get_membership(bc075_800_exp2, all_regions_BT075_800)
TLDR_idx_BT075_800 <- get_membership(bc075_800_TLDR, all_regions_BT075_800)
SLIC_idx_BT075_800 <- get_membership(bc075_800_SLIC, all_regions_BT075_800)

# Convert to character
makeVenn(list(
  "Exp1 BT075 W800" = as.character(exp1_idx_BT075_800),
  "Exp2 BT075 W800" = as.character(exp2_idx_BT075_800),
  "TLDR BT075 W800" = as.character(TLDR_idx_BT075_800),
  "SLIC BT075 W800" = as.character(SLIC_idx_BT075_800)
))






# Find overlaps stepwise
hits_1v2 <- findOverlaps(bc075_400_exp1, bc075_400_exp2, ignore.strand=TRUE)
gr_1v2 <- pintersect(
  bc075_400_exp1[queryHits(hits_1v2)],
  bc075_400_exp2[subjectHits(hits_1v2)]
)

hits_12v3 <- findOverlaps(gr_1v2, bc075_400_TLDR, ignore.strand=TRUE)
gr_12v3 <- pintersect(
  gr_1v2[queryHits(hits_12v3)],
  bc075_400_TLDR[subjectHits(hits_12v3)]
)

hits_123v4 <- findOverlaps(gr_12v3, bc075_400_SLIC, ignore.strand=TRUE)
gr_123v4 <- pintersect(
  gr_12v3[queryHits(hits_123v4)],
  bc075_400_SLIC[subjectHits(hits_123v4)]
)

export(gr_123v4, "common_all4_full.bed", format = "BED")

#Extend since we do not have TSS
gr_extended <- resize(gr_123v4, width = width(gr_123v4) + 100, fix = "center")

export(gr_extended, "common_all4_full_extended.bed", format = "BED")



