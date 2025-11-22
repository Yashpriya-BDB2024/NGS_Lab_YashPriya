library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(scales)

# Read your file
df <- read_tsv("subset_genera.bracken")

# Extract Genus from 'name' (first word)
df <- df %>%
  mutate(Genus = word(name, 1))

# Sum fraction_total_reads by Genus
df_genus <- df %>%
  group_by(Genus) %>%
  summarise(total_fraction = sum(fraction_total_reads)) %>%
  ungroup()

# Add SampleID (since you have only one sample)
df_genus <- df_genus %>%
  mutate(SampleID = "Sample1")

# 👉 Group low-abundance genera (<3%) into 'Other'
df_genus <- df_genus %>%
  mutate(Genus = ifelse(total_fraction < 0.03, "Other", Genus)) %>%
  group_by(Genus, SampleID) %>%
  summarise(total_fraction = sum(total_fraction)) %>%
  ungroup()

# Check grouped table
print(df_genus %>% arrange(desc(total_fraction)))

# ✅ Plot
p_all <- ggplot(df_genus, aes(x = SampleID, y = total_fraction, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Genus Relative Abundance (≥3% shown; others grouped)",
    x = "Sample ID",
    y = "Relative Abundance (%)"
  ) +
  scale_fill_manual(values = scales::hue_pal()(length(unique(df_genus$Genus)))) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

# Display plot
print(p_all)

# Save to file
ggsave("Genus_Relative_Abundance_Grouped.png", p_all, width = 10, height = 8, dpi = 300)
