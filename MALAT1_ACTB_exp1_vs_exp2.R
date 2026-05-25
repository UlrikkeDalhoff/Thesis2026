
library(ggplot2)

Malat1_exp1_PolyAplus <- read.table("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/PolyA+:PolyA-/MALAT1_real_PolyAplus_exp1.bed", sep = "\t", header = FALSE)
Malat1_exp1_PolyAminus <- read.table("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/PolyA+:PolyA-/MALAT1_real_PolyAminus_exp1.bed", sep = "\t", header = FALSE)

Malat1_exp1_PolyAplus$source <- "PolyA_plus_exp1"
Malat1_exp1_PolyAminus$source <- "PolyA_minus_exp1"

merged_MALAT1_Exp1 <- rbind(Malat1_exp1_PolyAminus, Malat1_exp1_PolyAplus) 

Malat1_exp2_PolyAplus <- read.table("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/PolyA+:PolyA-/MALAT1_real_PolyAplus_exp2.bed", sep = "\t", header = FALSE)
Malat1_exp2_PolyAminus <- read.table("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/PolyA+:PolyA-/MALAT1_real_PolyAminus_exp2.bed", sep = "\t", header = FALSE)

Malat1_exp2_PolyAplus$source <- "PolyA_plus_exp2"
Malat1_exp2_PolyAminus$source <- "PolyA_minus_exp2"

merged_MALAT1_Exp2 <- rbind(Malat1_exp2_PolyAminus, Malat1_exp2_PolyAplus) 

ACTB_exp1_PolyAplus <- read.table("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/PolyA+:PolyA-/ACTB_real_PolyAplus_exp1.bed", sep = "\t", header = FALSE)
ACTB_exp1_PolyAminus <- read.table("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/PolyA+:PolyA-/ACTB_real_PolyAminus_exp1.bed", sep = "\t", header = FALSE)

ACTB_exp1_PolyAplus$source <- "PolyA_plus_exp1"
ACTB_exp1_PolyAminus$source <- "PolyA_minus_exp1"

merged_ACTB_Exp1 <- rbind(ACTB_exp1_PolyAminus, ACTB_exp1_PolyAplus) 

ACTB_exp2_PolyAplus <- read.table("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/PolyA+:PolyA-/ACTB_real_PolyAplus_exp2.bed", sep = "\t", header = FALSE)
ACTB_exp2_PolyAminus <- read.table("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/PolyA+:PolyA-/ACTB_real_PolyAminus_exp2.bed", sep = "\t", header = FALSE)

ACTB_exp2_PolyAplus$source <- "PolyA_plus_exp2"
ACTB_exp2_PolyAminus$source <- "PolyA_minus_exp2"

merged_ACTB_Exp2 <- rbind(ACTB_exp2_PolyAminus, ACTB_exp2_PolyAplus) 



#visualize it 
# ggplot(merged_MALAT1_Exp1, aes(x = V11, fill = source)) +
#   geom_histogram(alpha = 0.5, position = "identity", bins = 50) +
#   theme_minimal() +
#   ggtitle("MALAT1 fractions of modifications in experiment 1") +
#   scale_x_continuous(breaks = seq(5,100, by = 10)) +
#   labs(x = "Fraction of modification", y = "Count of reads")
# 
# ggplot(merged_MALAT1_Exp2, aes(x = V11, fill = source)) +
#   geom_histogram(alpha = 0.5, position = "identity", bins = 50) +
#   theme_minimal() +
#   scale_x_continuous(breaks = seq(5,100, by = 10)) +
#   ggtitle("MALAT1 fractions of modifications in experiment 2") +
#   labs(x = "Fraction of modification", y = "Count of reads")
# 
# ggplot(merged_ACTB_Exp1, aes(x = V11, fill = source)) +
#   geom_histogram(alpha = 0.5, position = "identity", bins = 50) +
#   theme_minimal() +
#   ggtitle("ACTB fractions of modifications in experiment 1") +
#   scale_x_continuous(breaks = seq(5,100, by = 10)) +
#   labs(x = "Fraction of modification", y = "Count of reads")
# 
# ggplot(merged_ACTB_Exp2, aes(x = source, y = V11)) +
#   geom_col() +
#   theme_minimal() +
#   ggtitle("ACTB fractions of modifications in experiment 2") +
#   labs(x = "Fraction of modification", y = "Count of reads")


sites <- sort(unique(merged_ACTB_Exp2$V2))
breaks_to_show <- sites[seq(1, length(sites), by = 20)]

ggplot(merged_ACTB_Exp2, aes(x = factor(V2), y = V11, fill = source)) +
  geom_col(position = "dodge") +
  scale_x_discrete(breaks = breaks_to_show) + 
  labs(x = "Site", y = "Fraction modified", fill = "PolyA status") +
  ggtitle("ACTB fractions of modifications for each site in experiment 2") +
  theme_minimal() 
  
sites_2 <- sort(unique(merged_ACTB_Exp1$V2))
breaks_to_show_2 <- sites_2[seq(1, length(sites_2), by = 20)]

ggplot(merged_ACTB_Exp1, aes(x = factor(V2), y = V11, fill = source)) +
  geom_col(position = "dodge") +
  scale_x_discrete(breaks = breaks_to_show_2) +
  labs(x = "Site", y = "Fraction modified", fill = "PolyA status") +
  ggtitle("ACTB fractions of modifications for each site in experiment 1") +
  theme_minimal()

sites_3 <- sort(unique(merged_MALAT1_Exp2$V2))
breaks_to_show_3 <- sites_3[seq(1, length(sites_3), by = 20)]

ggplot(merged_MALAT1_Exp2, aes(x = factor(V2), y = V11, fill = source)) +
  geom_col(position = "dodge") +
  scale_x_discrete(breaks = breaks_to_show_3) +
  labs(x = "Site", y = "Fraction modified", fill = "PolyA status") +
  ggtitle("MALAT1 fractions of modifications for each site in experiment 2") +
  theme_minimal()

sites_4 <- sort(unique(merged_MALAT1_Exp1$V2))
breaks_to_show_4 <- sites_4[seq(1, length(sites_4), by = 20)]

ggplot(merged_MALAT1_Exp1, aes(x = factor(V2), y = V11, fill = source)) +
  geom_col(position = "dodge") +
  scale_x_discrete(breaks = breaks_to_show_4) +
  labs(x = "Site", y = "Fraction modified", fill = "PolyA status") +
  ggtitle("MALAT1 fractions of modifications for each site in experiment 1") +
  theme_minimal()

