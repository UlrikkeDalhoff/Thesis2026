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

#example use case - just give it lists of ids/names
makeVenn(list("list_1"=c("A","B","1_493","9_2022"),
              "list_2"=c("C","B","1_493","6_391","9_2022","1_292")))



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
bc095_100_exp1 <- clusterBidirectionally(
  ctss_exp1,
  balanceThreshold = 0.95,
  window = 101
)

bc095_200_exp1 <- clusterBidirectionally(
  ctss_exp1,
  balanceThreshold = 0.95,
  window = 201
)

bc095_400_exp1 <- clusterBidirectionally(
  ctss_exp1,
  balanceThreshold = 0.95,
  window = 401
)

bc095_600_exp1 <- clusterBidirectionally(
  ctss_exp1,
  balanceThreshold = 0.95,
  window = 601
)

bc095_800_exp1 <- clusterBidirectionally(
  ctss_exp1,
  balanceThreshold = 0.95,
  window = 801
)


#Now we can detect bidirectional clusteres (enhancers)
bc095_100_exp2 <- clusterBidirectionally(
  ctss_exp2,
  balanceThreshold = 0.95,
  window = 101
)

bc095_200_exp2 <- clusterBidirectionally(
  ctss_exp2,
  balanceThreshold = 0.95,
  window = 201
)

bc095_400_exp2 <- clusterBidirectionally(
  ctss_exp2,
  balanceThreshold = 0.95,
  window = 401
)

bc095_600_exp2 <- clusterBidirectionally(
  ctss_exp2,
  balanceThreshold = 0.95,
  window = 601
)

bc095_800_exp2 <- clusterBidirectionally(
  ctss_exp2,
  balanceThreshold = 0.95,
  window = 801
)

#Now we can detect bidirectional clusteres (enhancers)
bc095_100_TLDR <- clusterBidirectionally(
  ctss_TLDR,
  balanceThreshold = 0.95,
  window = 101
)

bc095_200_TLDR <- clusterBidirectionally(
  ctss_TLDR,
  balanceThreshold = 0.95,
  window = 201
)

bc095_400_TLDR <- clusterBidirectionally(
  ctss_TLDR,
  balanceThreshold = 0.95,
  window = 401
)

bc095_600_TLDR <- clusterBidirectionally(
  ctss_TLDR,
  balanceThreshold = 0.95,
  window = 601
)

bc095_800_TLDR <- clusterBidirectionally(
  ctss_TLDR,
  balanceThreshold = 0.95,
  window = 801
)


#Now we can detect bidirectional clusteres (enhancers)
bc095_100_SLIC <- clusterBidirectionally(
  ctss_SLIC,
  balanceThreshold = 0.95,
  window = 101
)

bc095_200_SLIC <- clusterBidirectionally(
  ctss_SLIC,
  balanceThreshold = 0.95,
  window = 201
)

bc095_400_SLIC <- clusterBidirectionally(
  ctss_SLIC,
  balanceThreshold = 0.95,
  window = 401
)

bc095_600_SLIC <- clusterBidirectionally(
  ctss_SLIC,
  balanceThreshold = 0.95,
  window = 601
)

bc095_800_SLIC <- clusterBidirectionally(
  ctss_SLIC,
  balanceThreshold = 0.95,
  window = 801
)


library(GenomicRanges)
library(eulerr)


# Combine all regions
all_regions <- reduce(c(
  bc095_100_exp1,
  bc095_100_exp2,
  bc095_100_TLDR,
  bc095_100_SLIC
))

# Function to check overlap membership
get_membership <- function(query, subject) {
  hits <- findOverlaps(subject, query, ignore.strand = TRUE)
  unique(queryHits(hits))
}

# Membership indices in common universe 
exp1_idx <- get_membership(bc095_100_exp1, all_regions)
exp2_idx <- get_membership(bc095_100_exp2, all_regions)
TLDR_idx <- get_membership(bc095_100_TLDR, all_regions)
SLIC_idx <- get_membership(bc095_100_SLIC, all_regions)

# Convert to character
makeVenn(list(
  "Exp1 BT095 W100" = as.character(exp1_idx),
  "Exp2 BT095 W100" = as.character(exp2_idx),
  "TLDR BT095 W100" = as.character(TLDR_idx),
  "SLIC BT095 W100" = as.character(SLIC_idx)
))

#for BT095 window 200
# Combine all regions
all_regions_BT095_200 <- reduce(c(
  bc095_200_exp1,
  bc095_200_exp2,
  bc095_200_TLDR,
  bc095_200_SLIC
))


# Membership indices in common universe 
exp1_idx_BT095_200 <- get_membership(bc095_200_exp1, all_regions_BT095_200)
exp2_idx_BT095_200 <- get_membership(bc095_200_exp2, all_regions_BT095_200)
TLDR_idx_BT095_200 <- get_membership(bc095_200_TLDR, all_regions_BT095_200)
SLIC_idx_BT095_200 <- get_membership(bc095_200_SLIC, all_regions_BT095_200)

# Convert to character
makeVenn(list(
  "Exp1 BT095 W200" = as.character(exp1_idx_BT095_200),
  "Exp2 BT095 W200" = as.character(exp2_idx_BT095_200),
  "TLDR BT095 W200" = as.character(TLDR_idx_BT095_200),
  "SLIC BT095 W200" = as.character(SLIC_idx_BT095_200)
))

#for BT095 window 400
# Combine all regions
all_regions_BT095_400 <- reduce(c(
  bc095_400_exp1,
  bc095_400_exp2,
  bc095_400_TLDR,
  bc095_400_SLIC
))


# Membership indices in common universe 
exp1_idx_BT095_400 <- get_membership(bc095_400_exp1, all_regions_BT095_400)
exp2_idx_BT095_400 <- get_membership(bc095_400_exp2, all_regions_BT095_400)
TLDR_idx_BT095_400 <- get_membership(bc095_400_TLDR, all_regions_BT095_400)
SLIC_idx_BT095_400 <- get_membership(bc095_400_SLIC, all_regions_BT095_400)

# Convert to character
makeVenn(list(
  "Exp1 BT095 W400" = as.character(exp1_idx_BT095_400),
  "Exp2 BT095 W400" = as.character(exp2_idx_BT095_400),
  "TLDR BT095 W400" = as.character(TLDR_idx_BT095_400),
  "SLIC BT095 W400" = as.character(SLIC_idx_BT095_400)
))




#for BT095 window 600
# Combine all regions
all_regions_BT095_600 <- reduce(c(
  bc095_600_exp1,
  bc095_600_exp2,
  bc095_600_TLDR,
  bc095_600_SLIC
))


# Membership indices in common universe 
exp1_idx_BT095_600 <- get_membership(bc095_600_exp1, all_regions_BT095_600)
exp2_idx_BT095_600 <- get_membership(bc095_600_exp2, all_regions_BT095_600)
TLDR_idx_BT095_600 <- get_membership(bc095_600_TLDR, all_regions_BT095_600)
SLIC_idx_BT095_600 <- get_membership(bc095_600_SLIC, all_regions_BT095_600)

# Convert to character
makeVenn(list(
  "Exp1 BT095 W600" = as.character(exp1_idx_BT095_600),
  "Exp2 BT095 W600" = as.character(exp2_idx_BT095_600),
  "TLDR BT095 W600" = as.character(TLDR_idx_BT095_600),
  "SLIC BT095 W600" = as.character(SLIC_idx_BT095_600)
))



#for BT095 window 800
# Combine all regions
all_regions_BT095_800 <- reduce(c(
  bc095_800_exp1,
  bc095_800_exp2,
  bc095_800_TLDR,
  bc095_800_SLIC
))


# Membership indices in common universe 
exp1_idx_BT095_800 <- get_membership(bc095_800_exp1, all_regions_BT095_800)
exp2_idx_BT095_800 <- get_membership(bc095_800_exp2, all_regions_BT095_800)
TLDR_idx_BT095_800 <- get_membership(bc095_800_TLDR, all_regions_BT095_800)
SLIC_idx_BT095_800 <- get_membership(bc095_800_SLIC, all_regions_BT095_800)

# Convert to character
makeVenn(list(
  "Exp1 BT095 W800" = as.character(exp1_idx_BT095_800),
  "Exp2 BT095 W800" = as.character(exp2_idx_BT095_800),
  "TLDR BT095 W800" = as.character(TLDR_idx_BT095_800),
  "SLIC BT095 W800" = as.character(SLIC_idx_BT095_800)
))