library(tidyverse)
header <- c("chr","start","end","count")
df19 <- read_tsv("/Users/cmdb/qb25-answers/week1/hg19-kc-count.bed", col_names = header)
df16 <- read_tsv("/Users/cmdb/qb25-answers/week1/hg16-kc-count.bed", col_names = header)
df <- bind_rows(list(hg19 = df19, hg16 = df16), .id = "assembly")

ggplot(df, aes(x = start, y = count, color = assembly, group = interaction (assembly, chr))) +
  geom_point(size = 0.3) + 
  geom_line() +
  scale_y_log10() +
  labs(title = "Gene density",
       x = "Start",
       y = "Count") +
  facet_wrap(~ chr, scales = "free")
ggsave("exercise2.png", width = 12, height = 10, dpi = 300)

