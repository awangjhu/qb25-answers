library(tidyverse)

#AF Histogram
af_lines <- read_lines("AF.txt")
af_lines <- af_lines[af_lines != "AF" & af_lines != ""]
af <- data.frame(AF = as.numeric(af_lines))
af <- subset(af, !is.na(AF))

p_af <- ggplot(af, aes(x = AF)) +
  geom_histogram(bins = 11) +
  labs(
    title = "Allele Frequency Spectrum (BY×RM segregants)",
    x = "Allele frequency (ALT across samples)",
    y = "Variant count"
  ) +
  theme_classic()

ggsave("AF_hist.png", p_af, width = 8, height = 5, dpi = 150)

#Interpretation of the figure:
#The histogram is showing the distribution of allele frequencies from the samples.
#It seems that more variants are clustered around the middle frequencies (0.25-0.75) rather than near the extremes,
#likely because the allele segregates evenly. This distribution is called the allele frequency spectrum and the shape is 
#supposed to resemble a binomial distribution.

#DP Histogram
lines <- readr::read_lines("DP.txt")
vals  <- suppressWarnings(as.numeric(lines))
dp    <- data.frame(DP = vals)
dp    <- subset(dp, !is.na(DP) & is.finite(DP) & DP >= 0)

p <- ggplot(dp, aes(x = DP)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left") +
  coord_cartesian(xlim = c(0, 20)) +
  scale_x_continuous(breaks = 0:20) +
  labs(
    title = "Read Depth Distribution Across All Sample–Variant Genotypes",
    x = "Read depth (DP)", y = "Count"
  ) +
  theme_classic()

ggsave("DP_hist.png", p, width = 8, height = 5, dpi = 150)

#Interpretation of the figure
#This histogram is showing how many variant fall into each read depth across the samples. Based on the histogram,
#it appears that most sites have a read depth between 2-6. Yes, this is expected because the data comes from illumina short-read sequencing.
#I believe the distribution is a poisson distribution.


