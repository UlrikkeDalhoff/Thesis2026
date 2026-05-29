library(CAGEfightR)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)

#We create a CTSS object and get a count at each 5' position
ctss <- quantifyCTSSs(
  plusStrand = c(sample1 = "/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancer_supportdata/TLDR_plus_clean2.bedGraph"),
  minusStrand = c(sample1 = "/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/enhancer_supportdata/TLDR_minus_clean2.bedGraph"),
  genome = seqinfo(BSgenome.Hsapiens.UCSC.hg38)
)


#And then we normalize
ctss <- calcTPM(ctss)

#Combine samples 
ctss <- calcPooled(ctss)

#Now we can detect bidirectional clusteres (enhancers)
bc095_100 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.95,
  window = 101
)

bc095_200 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.95,
  window = 201
)

bc095_400 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.95,
  window = 401
)

bc095_600 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.95,
  window = 601
)

bc095_800 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.95,
  window = 801
)


#bt = balancethreshold. 
window <- c("100", "200","400", "600", "800")
count <- c(length(bc095_100$thick), length(bc095_200$thick), length(bc095_400$thick), length(bc095_600$thick), length(bc095_800$thick))

BT095 <- data.frame(window, count)



ggplot(BT095, aes(x = window, y = count, fill = window)) +
  geom_col(show.legend = FALSE) +
  labs(
    title = "TLDR data with Balancethreshold = 0.95",
    x = "window size (bp)",
    y = "Count of bidirectional clusters"
  ) +
  theme_minimal()












#Now we can detect bidirectional clusteres (enhancers)
bc050_100 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.50,
  window = 101
)

bc050_200 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.50,
  window = 201
)

bc050_400 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.50,
  window = 401
)

bc050_600 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.50,
  window = 601
)

bc050_800 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.50,
  window = 801
)


#bt = balancethreshold. 
window_50 <- c("100", "200","400", "600", "800")
count_50 <- c(length(bc050_100$thick), length(bc050_200$thick), length(bc050_400$thick), length(bc050_600$thick), length(bc050_800$thick))

BT050 <- data.frame(window_50, count_50)



ggplot(BT050, aes(x = window_50, y = count_50, fill = window_50)) +
  geom_col(show.legend = FALSE) +
  labs(
    title = "TLDR data with Balancethreshold = 0.50",
    x = "window size (bp)",
    y = "Count of bidirectional clusters"
  ) +
  theme_minimal()









#Now we can detect bidirectional clusteres (enhancers)
bc075_100 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.75,
  window = 101
)

bc075_200 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.75,
  window = 201
)

bc075_400 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.75,
  window = 401
)

bc075_600 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.75,
  window = 601
)

bc075_800 <- clusterBidirectionally(
  ctss,
  balanceThreshold = 0.75,
  window = 801
)


#bt = balancethreshold. 
window_75 <- c("100", "200","400", "600", "800")
count_75 <- c(length(bc075_100$thick), length(bc075_200$thick), length(bc075_400$thick), length(bc075_600$thick), length(bc075_800$thick))

BT075 <- data.frame(window_75, count_75)



ggplot(BT075, aes(x = window_75, y = count_75, fill = window_75)) +
  geom_col(show.legend = FALSE) +
  labs(
    title = "TLDR data with Balancethreshold = 0.75",
    x = "window size (bp)",
    y = "Count of bidirectional clusters"
  ) +
  theme_minimal()

