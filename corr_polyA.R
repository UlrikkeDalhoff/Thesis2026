polyAplus_bed <- read.table("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/merged_data_polyA/m6Arealresult2_polyAplus.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
colnames(polyAplus_bed)[1] <- "Chromosome"
colnames(polyAplus_bed)[2] <- "Start"
colnames(polyAplus_bed)[3] <- "End"
colnames(polyAplus_bed)[4] <- "Mod"
colnames(polyAplus_bed)[5] <- "Reads"
colnames(polyAplus_bed)[6] <- "Strand"
colnames(polyAplus_bed)[11] <- "Frac mod"
colnames(polyAplus_bed)[12] <- "N mod"


polyAminus_bed <- read.table("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/merged_data_polyA/m6Arealresult2_polyAminus.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
colnames(polyAminus_bed)[1] <- "Chromosome"
colnames(polyAminus_bed)[2] <- "Start"
colnames(polyAminus_bed)[3] <- "End"
colnames(polyAminus_bed)[4] <- "Mod"
colnames(polyAminus_bed)[5] <- "Reads"
colnames(polyAminus_bed)[6] <- "Strand"
colnames(polyAminus_bed)[11] <- "Frac mod"
colnames(polyAminus_bed)[12] <- "N mod"

library(ggplot2)
library(tidyverse)

polyAplus_bed$Condition <- "polyA_plus"
polyAminus_bed$Condition <- "polyA_minus"

combined <- rbind(polyAplus_bed, polyAminus_bed)

combined$Reads <- as.numeric(combined$Reads)
combined$`Frac mod` <- as.numeric(combined$`Frac mod`)


#ggplot(combined, aes(x = Reads, y = `Frac mod`, color = Condition)) +
#  geom_point(alpha = 0.5) +
#  theme_classic() +
#  scale_x_log10() +
#  labs(x = "Reads", y = "Fraction modified")


#ggplot(combined, aes(x = Reads, y = `Frac mod`)) +
#  geom_point(alpha = 0.5) +
#  facet_wrap(~Condition) +
#  theme_classic() +
#  scale_x_log10() +
#  labs(x = "Reads", y = "Fraction modified")

#ggplot(combined, aes(x = Reads, fill = Condition)) +
#  geom_density(alpha = 0.4) +
#  scale_x_log10() +
#  theme_classic()

#ggplot(combined, aes(x = `Frac mod`, fill = Condition)) +
#  geom_density(alpha = 0.4) +
#  scale_x_log10() +
#  theme_classic()


#For correlationn
merged <- merge(
  polyAplus_bed,
  polyAminus_bed,
  by = c("Chromosome", "Start", "End"),
  suffixes = c("_plus", "_minus")
)

ggplot(merged, aes(x = `Frac mod_plus`, y = `Frac mod_minus`)) +
  geom_point(alpha = 0.5) +
  geom_abline() +
  theme_classic() +
  labs(x = "Frac mod polyA+", y = "Frac mod polyA-")

ggplot(merged, aes(x = `N mod_plus`, y = `N mod_minus`)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(x = "N mod polyA+", y = "N mod polyA-")


r_val <- cor(merged$`Frac mod_plus`, merged$`Frac mod_minus`, 
             use = "complete.obs", method = "pearson")

ggplot(merged, aes(x = `Frac mod_plus`, y = `Frac mod_minus`)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_classic() +
  ggtitle("Scatter plot of mods in PolyA- vs PolyA+") +
  labs(x = "Frac mod polyA+", y = "Frac mod polyA-") +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("Pearson's r = ", round(r_val, 3)),
           hjust = -0.1, vjust = 1.5, size = 4)




ggplot(merged, aes(x = `Frac mod_plus`, y = `Frac mod_minus`)) +
  geom_point(alpha = 0.5) +
  geom_abline(stat = "identity", col = "red") +
  theme_classic() +
  ggtitle("Scatter plot of mods in PolyA- vs PolyA+") +
  labs(x = "Frac mod polyA+", y = "Frac mod polyA-")



merged <- merged %>%
  mutate(row_id = row_number())

library(plotly)

p <- ggplot(merged, aes(x = `Frac mod_plus`, y = `Frac mod_minus`,
                        text = row_id)) +  # shows on hover
  geom_point(alpha = 0.5) +
  geom_abline(col = "red") +
  theme_classic() +
  labs(x = "Frac mod polyA+", y = "Frac mod polyA-")

ggplotly(p, tooltip = "text")


#2569, 3495, 790, 4086, 4124 and 4933 

Outliers <- merged %>%
  filter(row_id == "2569" | row_id == "3495" | row_id == "790" | row_id == "4086" | row_id == "4124" | row_id == "4933")


write.table(
  Outliers,
  file = "outliers_Real_poly.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)



