library(tidyverse)
library(scales)

#visualizing one sample (3.3)
gt <- read_tsv("/Users/cmdb/qb25-answers/week3/gt_long.txt")

# filter to one sample and one chromosome
df <- subset(gt, sample == "A01_62" & chrom == "chrII")

ggplot(df, aes(x = as.numeric(pos), y = factor(genotype))) +
  geom_point(alpha = 0.6) +
  labs(
    title = "Ancestry of sample A01_62 on chrII",
    x = "Position on chrII",
    y = "Genotype (0=ref, 1=alt)"
  ) +
  theme_classic()

#Interpretation:
#In sample A01_62, for chromosome II, there tracts of consecutive or alternative genotypes, shown as 0 or 1 respectively. These mark crossover points between the BY lab
#strain and the RM wine strain ancestory. This is what I expect in mosaic inheritance.

df1 <- subset(gt, sample == "A01_62")


gt <- read_tsv("gt_long.txt", show_col_types = FALSE) |>
  mutate(pos = as.numeric(pos),
         Mb  = pos / 1e6,
         genotype = factor(genotype))

p_all <- ggplot(gt, aes(x = Mb, y = sample, color = genotype)) +
  geom_point(alpha = 0.6, size = 0.4) +
  facet_grid(~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous("Position (Mb)",
                     breaks = breaks_pretty(n = 3),
                     labels = label_number(accuracy = 0.1)) +
  labs(title = "Ancestry across all samples", y = "Sample", color = "Genotype") +
  theme_classic(base_size = 11) +
  theme(panel.spacing.x = unit(0.6, "lines"),
        axis.text.x = element_text(angle = 0, vjust = 1))

setwd("/Users/cmdb/qb25-answers/week3")
ggsave("All_samples_ancestry_facets.png", p_all, width = 12, height = 6, dpi = 200)
