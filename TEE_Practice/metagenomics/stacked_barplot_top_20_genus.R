library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(scales)  # for hue_pal() and percent labels

# 1️⃣ Read and preprocess
df <- read_tsv("subset_genera.bracken")

# Extract Genus from 'name'
df <- df %>%
  mutate(Genus = word(name, 1))

# Sum fraction_total_reads by Genus
df_genus <- df %>%
  group_by(Genus) %>%
  summarise(total_fraction = sum(fraction_total_reads)) %>%
  ungroup()

# Add a single sample ID
df_genus <- df_genus %>%
  mutate(SampleID = "Sample1")

# 2️⃣ Identify Top 20 Genera
topN <- 20

top_genera <- df_genus %>%
  top_n(topN, total_fraction) %>%
  pull(Genus)

# Group others as "Other"
df_top20 <- df_genus %>%
  mutate(Genus_grouped = ifelse(Genus %in% top_genera, Genus, "Other")) %>%
  group_by(SampleID, Genus_grouped) %>%
  summarise(total_fraction = sum(total_fraction), .groups = "drop")

# 3️⃣ Plot Top 20 + Other
p_top20 <- ggplot(df_top20, aes(x = SampleID, y = total_fraction, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = paste("Genus Relative Abundance (Top", topN, "+ Other)"),
    x = "Sample ID",
    y = "Relative Abundance (%)"
  ) +
  scale_fill_manual(values = scales::hue_pal()(length(unique(df_top20$Genus_grouped)))) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

# 4️⃣ Display and Save
print(p_top20)

ggsave("Genus_Relative_Abundance_Top20.png", p_top20, width = 10, height = 8, dpi = 300)
