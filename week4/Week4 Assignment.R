#Exercise 1:
library(tidyverse)
library(broom)
#1.1:
dnm <- read_csv("/Users/cmdb/qb25-answers/week4/aau1043_dnm.csv")
glimpse(dnm)
head(dnm)

#1.2: 
dnm_counts <- dnm %>%
  filter(!is.na(Phase_combined)) %>%
  group_by(Proband_id, Phase_combined) %>%
  summarise(dnm_count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = Phase_combined,
    values_from = dnm_count,
    values_fill = 0
  ) %>%
  rename(
    maternal_dnm = mother,
    paternal_dnm = father
  )
head(dnm_counts)
#1.3/1.4: 
ages <- read_csv("/Users/cmdb/qb25-answers/week4/aau1043_parental_age.csv") %>%
  rename(
    paternal_age = Father_age,
    maternal_age = Mother_age
  )
glimpse(ages)
#Merge count w/ age
merged <- dnm_counts %>%
  inner_join(ages, by = "Proband_id")

glimpse(merged)
head(merged)

#Exercise 2:
#2.1
p_maternal <- merged %>%
  ggplot(aes(x = maternal_age, y = maternal_dnm)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  labs(
    x = "Maternal Age",
    y = "Maternal DNM Count",
    title = "Maternal DNM Count vs. Maternal Age"
  ) +
  theme_minimal()
setwd("/Users/cmdb/qb25-answers/week4")
ggsave("ex2_a.png", p_maternal, width = 6, height = 4, dpi = 300)

p_paternal <- merged %>%
  ggplot(aes(x = paternal_age, y = paternal_dnm)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  labs(
    x = "Paternal Age",
    y = "Paternal DNM Count",
    title = "Paternal DNM Count vs. Paternal Age"
  ) +
  theme_minimal()
setwd("/Users/cmdb/qb25-answers/week4")
ggsave("ex2_b.png", p_paternal, width = 6, height = 4, dpi = 300)
#2.2
model_maternal <- lm(maternal_dnm ~ maternal_age, data = merged)
maternal_results <- broom::tidy(model_maternal)
maternal_results
#2.3
model_paternal <- lm(paternal_dnm ~ paternal_age, data = merged)
paternal_results <- broom::tidy(model_paternal)
paternal_results
#2.4
new_data <- tibble(paternal_age = 50.5)
prediction <- predict(model_paternal, new_data, interval = "confidence")
prediction
#2.6:
#make easier for ploting
dnm_long <- merged %>%
  pivot_longer(
    cols = c(maternal_dnm, paternal_dnm),
    names_to = "origin",
    values_to = "dnm_count"
  )
p_hist <- dnm_long %>%
  ggplot(aes(x = dnm_count, fill = origin)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(
    x = "DNM Count",
    y = "Frequency",
    title = "Distribution of Maternal vs. Paternal DNMs"
  ) +
  scale_fill_manual(values = c("maternal_dnm" = "tomato", "paternal_dnm" = "steelblue")) +
  theme_minimal()
setwd("/Users/cmdb/qb25-answers/week4")
ggsave("ex2_c.png", p_hist, width = 6, height = 4, dpi = 300)
#2.6:
#Paired t-test
t_test_result <- t.test(
  merged$maternal_dnm,
  merged$paternal_dnm,
  paired = TRUE
)

t_test_result

#Exercise 3
#3.1: I picked pokemon
pokemon_df <- read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/main/data/2025/2025-04-01/pokemon_df.csv')
glimpse(pokemon_df)
#3.2: 
#First, I wawnt to test if speed and attack are inversely related in Pokemon
model_speed_attack <- lm(speed ~ attack, data = pokemon_df)
broom::tidy(model_speed_attack)
p_speed_attack <- pokemon_df %>%
  ggplot(aes(x = attack, y = speed)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  labs(
    title = "Speed vs. Attack Tradeoff in Pokémon",
    x = "Attack",
    y = "Speed"
  ) +
  theme_minimal()
setwd("/Users/cmdb/qb25-answers/week4")
ggsave("ex3_speed_attack.png", p_speed_attack, width = 6, height = 4, dpi = 300)

#Second I want to test if total base stats correlate from weight or type
pokemon_df <- pokemon_df %>%
  mutate(total_stats = hp + attack + defense + special_attack + special_defense + speed)
#total stats by log weight
model_stats_weight <- lm(total_stats ~ log10(weight), data = pokemon_df)
broom::tidy(model_stats_weight)
#Total stats by primary type
model_stats_type <- lm(total_stats ~ type_1, data = pokemon_df)
broom::tidy(model_stats_type)

p_stats_weight <- pokemon_df %>%
  mutate(total_stats = hp + attack + defense + special_attack + special_defense + speed) %>%
  ggplot(aes(x = log10(weight), y = total_stats)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  labs(
    title = "Total Base Stats vs. Pokémon Weight",
    x = "Weight (log10 scale)",
    y = "Total Base Stats"
  ) +
  theme_minimal()
setwd("/Users/cmdb/qb25-answers/week4")
ggsave("ex3_totalstats_weight.png", p_stats_weight, width = 6, height = 4, dpi = 300)

#3.3: There is a relationship between Pokemon weight and total base stats. As in heavier pokemon tend to have higher base stats
model_weight_stats <- lm(total_stats ~ log10(weight), data = pokemon_df)
broom::tidy(model_weight_stats)


