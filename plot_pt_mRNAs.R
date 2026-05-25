
library(ggplot2)
library(dplyr)

MALAT1_minus <- read.delim("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/mRNA_polydata/PolyAminus_MALAT1.tsv", header = FALSE, sep = "\t")
MALAT1_minus$V25 <- as.numeric(sub(".*:", "", MALAT1_minus$V25))
MALAT1_minus <- MALAT1_minus %>% 
  filter(!is.na(V25))

MALAT1_plus <- read.delim("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/mRNA_polydata/PolyAplus_MALAT1.tsv", header = FALSE, sep = "\t")
MALAT1_plus$V25 <- as.numeric(sub(".*:", "", MALAT1_plus$V25))
MALAT1_plus <- MALAT1_plus %>% 
  filter(!is.na(V25))

ggplot(MALAT1_minus, aes(x = V25)) +
  geom_histogram(color = "lightblue", fill = "orange") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,15)) +
  labs(
    title = "Counts of Poly-A tail length in MALAT1 in Poly-A-",
    x = "Number of Nucleotides",
    y = "Counts")

ggplot(MALAT1_plus, aes(x = V25)) +
  geom_histogram(color = "orange", fill = "lightblue") +
  theme_classic() +
  scale_x_continuous(breaks = seq(11,400, by = 30)) +
  labs(
    title = "Counts of Poly-A tail length in MALAT1 in Poly-A+",
    x = "Number of Nucleotides",
    y = "Counts")


HIST1H4C_minus <- read.delim("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/mRNA_polydata/PolyAminus_HIST1H4C.tsv", header = FALSE, sep = "\t")
HIST1H4C_minus$V25 <- as.numeric(sub(".*:", "", HIST1H4C_minus$V25))
HIST1H4C_minus <- HIST1H4C_minus %>% 
  filter(!is.na(V25))

HIST1H4C_plus <- read.delim("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/mRNA_polydata/PolyAplus_HIST1H4C.tsv", header = FALSE, sep = "\t")
HIST1H4C_plus$V25 <- as.numeric(sub(".*:", "", HIST1H4C_plus$V25))
HIST1H4C_plus <- HIST1H4C_plus %>% 
  filter(!is.na(V25))

ggplot(HIST1H4C_minus, aes(x = V25)) +
  geom_histogram(color = "lightblue", fill = "orange") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,15)) +
  labs(
    title = "Density of Poly-A tail length in HIST1H4C in Poly-A-",
    x = "Number of Nucleotides",
    y = "Density")

ggplot(HIST1H4C_plus, aes(x = V25)) +
  geom_histogram(color = "orange", fill = "lightblue") +
  theme_classic() +
  scale_x_continuous(breaks = seq(11,400, by = 30)) +
  labs(
    title = "Density of Poly-A tail length in HIST1H4C in Poly-A+",
    x = "Number of Nucleotides",
    y = "Density")


B2M_minus <- read.delim("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/mRNA_polydata/PolyAminus_B2M.tsv", header = FALSE, sep = "\t")
B2M_minus$V25 <- as.numeric(sub(".*:", "", B2M_minus$V25))
B2M_minus <- B2M_minus %>% 
  filter(!is.na(V25))

B2M_plus <- read.delim("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/mRNA_polydata/PolyAplus_B2M.tsv", header = FALSE, sep = "\t")
B2M_plus$V25 <- as.numeric(sub(".*:", "", B2M_plus$V25))
B2M_plus <- B2M_plus %>% 
  filter(!is.na(V25))

ggplot(B2M_minus, aes(x = V25)) +
  geom_histogram(color = "lightblue", fill = "orange") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,15)) +
  labs(
    title = "Density of Poly-A tail length in B2M in Poly-A-",
    x = "Number of Nucleotides",
    y = "Density")

ggplot(B2M_plus, aes(x = V25)) +
  geom_histogram(color = "orange", fill = "lightblue") +
  theme_classic() +
  scale_x_continuous(breaks = seq(11,400, by = 30)) +
  labs(
    title = "Density of Poly-A tail length in B2M in Poly-A+",
    x = "Number of Nucleotides",
    y = "Density")

ACTB_plus <- read.delim("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/mRNA_polydata/PolyAplus_ACTB.tsv", header = FALSE, sep = "\t")
ACTB_plus$V25 <- as.numeric(sub(".*:", "", ACTB_plus$V25))
ACTB_plus <- ACTB_plus %>% 
  filter(!is.na(V25))

ACTB_minus <- read.delim("/Users/ulrikke/Desktop/Bioinformatics/Year2/Thesis/Data1/mRNA_polydata/PolyAminus_ACTB.tsv", header = FALSE, sep = "\t")
ACTB_minus$V25 <- as.numeric(sub(".*:", "", ACTB_minus$V25))
ACTB_minus <- ACTB_minus %>% 
  filter(!is.na(V25))

ggplot(ACTB_minus, aes(x = V25)) +
  geom_histogram(aes(y = after_stat(count)),
                 color = "orange", fill = "lightblue") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,15)) +
  labs(
    title = "Counts of Poly-A tail length in ACTB in Poly-A-",
    x = "Number of Nucleotides",
    y = "Counts")

ggplot(ACTB_plus, aes(x = V25)) +
  geom_histogram(aes(y = after_stat(count)),
                 color = "orange", fill = "lightblue") +
  theme_classic() +
  scale_x_continuous(breaks = seq(11,400, by = 30)) +
  labs(
    title = "Counts of Poly-A tail length in ACTB in Poly-A+",
    x = "Number of Nucleotides",
    y = "Counts")

