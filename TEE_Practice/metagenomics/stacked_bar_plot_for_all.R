library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# Read your file
df <- read_tsv("subset_species.bracken")

# Extract Species from 'name' (both first and second word)
df <- df %>%
  mutate(Species = word(name, 1, 2))

# Sum fraction_total_reads by Species
df_species <- df %>%
  group_by(Species) %>%
  summarise(total_fraction = sum(fraction_total_reads)) %>%
  ungroup()

# Add SampleID (since you have only one sample)
df_species <- df_species %>%
  mutate(SampleID = "Sample1")


### for all sample
library(ggplot2)
library(scales)   # for hue_pal()

p_all <- ggplot(df_species, aes(x = SampleID, y = total_fraction, fill = Species)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Species Relative Abundance (All Species)",
    x = "Sample ID",
    y = "Relative Abundance (%)"
  ) +
  scale_fill_manual(values = scales::hue_pal()(length(unique(df_species$Species)))) +
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
ggsave("Species_Relative_Abundance_All.png", p_all, width = 15, height = 8, dpi = 300)
