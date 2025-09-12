header <- c("chr", "start", "end", "count")
df_kc <- read_tsv("/Users/cmdb/qb25-answers/week1/hg19-kc-count.bed", col_names = header)
head(df_kc)

ggplot(df_kc, aes(x = start, y = count)) +
    geom_point(size = 0.3) + 
  geom_line() +
  scale_y_log10() +
  labs(title = "Gene density",
       x = "Start",
       y = "Count") +
  facet_wrap(~ chr, scales = "free")
ggsave("exercise1.png")