library(ggplot2)

# Read the bracken output
data <- read.table("subset_genera.bracken", header = TRUE, sep = "\t")

# Calculate richness and abundance
richness <- length(unique(data$name))
abundance <- sum(data$new_est_reads)

cat("Richness:", richness, "\n")
cat("Total Abundance (sum of reads):", abundance, "\n")

# Create a dummy dataframe for richness
richness_df <- data.frame(
  Metric = "Richness",
  Value = richness,
  x = "Overall"
)

# Combined plot with dual axes
pdf("diversity.pdf", width = 8, height = 6)

ggplot(data, aes(x = reorder(name, -new_est_reads))) +
  # Bar for abundance
  geom_bar(aes(y = new_est_reads), stat = "identity", fill = "steelblue") +
  
  # Line for richness (secondary axis)
  geom_line(
    data = richness_df,
    aes(x = x, y = Value * (max(data$new_est_reads) / richness)),
    color = "darkred",
    size = 1.2,
    group = 1
  ) +
  geom_point(
    data = richness_df,
    aes(x = x, y = Value * (max(data$new_est_reads) / richness)),
    color = "darkred",
    size = 3
  ) +
  
  scale_y_continuous(
    name = "Abundance (Reads)",
    sec.axis = sec_axis(~ . / (max(data$new_est_reads) / richness),
                        name = "Richness")
  ) +
  
  labs(
    x = "Taxa",
    title = "Species Abundance and Richness"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred")
  )

dev.off()

