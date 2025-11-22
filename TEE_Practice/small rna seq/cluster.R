# -----------------------------------------------------------------------------
# 1. Load All Necessary Libraries
# -----------------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)
library(ggplot2)
library(pathview) # For visualizing specific pathways

# -----------------------------------------------------------------------------
# 2. Define File Paths and Parameters
# -----------------------------------------------------------------------------

# Input file from your analysis (gene symbols and log2FC)
input_file <- "Human_DEG.txt"

# Output directory to save the results
output_dir <- "clusterprofiler_results"

# Minimum q-value for significant enrichment
qvalue_cutoff <- 0.05

# -----------------------------------------------------------------------------
# 3. Read Input Data and Prepare for Analysis
# -----------------------------------------------------------------------------

print("Step 1: Reading input data and converting to Entrez IDs...")
# Read the two-column TSV file
my_data <- read_tsv(input_file)

# Convert gene symbols to Entrez IDs
entrez_data <- bitr(
    my_data$gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)

# Join original data with Entrez IDs and remove unmapped genes
my_data_with_entrez <- my_data %>%
    left_join(entrez_data, by = c("gene" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID))

# Create a named vector of log2FC values
log2FC_vector <- my_data_with_entrez$log2FC
names(log2FC_vector) <- my_data_with_entrez$ENTREZID

# -----------------------------------------------------------------------------
# 4. Perform Enrichment Analysis (GO and KEGG)
# -----------------------------------------------------------------------------

# --- GO Enrichment ---
print("Step 2: Performing GO enrichment analysis...")
go_results <- enrichGO(
    gene = my_data_with_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = qvalue_cutoff,
    readable = TRUE # Get gene names in the result
)

# --- KEGG Enrichment ---
print("Step 3: Performing KEGG pathway enrichment...")
kegg_results <- enrichKEGG(
    gene = my_data_with_entrez$ENTREZID,
    organism = "hsa",  # "hsa" is the KEGG code for Homo sapiens
    pvalueCutoff = 0.05,
    qvalueCutoff = qvalue_cutoff
)

# -----------------------------------------------------------------------------
# 5. Save All Results (Data and Plots)
# -----------------------------------------------------------------------------

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

# --- Save GO data frame ---
print("Saving GO enrichment data table...")
write_tsv(
    as.data.frame(go_results),
    file.path(output_dir, "go_enrichment_results.tsv")
)

# --- Save KEGG data frame ---
print("Saving KEGG enrichment data table...")
write_tsv(
    as.data.frame(kegg_results),
    file.path(output_dir, "kegg_enrichment_results.tsv")
)


# --- Save GO Plots ---
print("Saving GO dot plot...")
p_dotplot <- dotplot(go_results, showCategory = 10, x = "GeneRatio")
ggsave(
    filename = "go_dotplot.png",
    plot = p_dotplot,
    path = output_dir,
    width = 8,
    height = 6,
    units = "in"
)

print("Saving GO gene-concept network plot with log2FC...")
top_go_term <- go_results@result$Description[1]
plot_title <- paste0("Gene-Concept Network Plot\nTop Pathway: ", top_go_term)

p_cnetplot <- cnetplot(
    go_results,
    showCategory = 5,
    foldChange = log2FC_vector,
    circular = FALSE
)

p_cnetplot_final <- p_cnetplot +
  ggplot2::theme_bw() +
  ggplot2::scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "log2FC"
  ) +
  ggplot2::labs(title = plot_title)

ggsave(
    filename = "go_cnetplot_log2FC.png",
    plot = p_cnetplot_final,
    path = output_dir,
    width = 10,
    height = 8,
    units = "in"
)

# --- Save KEGG Bar Plot ---
print("Saving KEGG bar plot...")
p_kegg_barplot <- barplot(kegg_results, showCategory = 10)
ggsave(
    filename = "kegg_barplot.png",
    plot = p_kegg_barplot,
    path = output_dir,
    width = 8,
    height = 6,
    units = "in"
)

# --- Save a Specific KEGG Pathway Diagram with log2FC ---
# You need to manually look up the KEGG ID for the pathway you want
target_kegg_id <- "hsa04110" # Example: Cell cycle pathway ID

print(paste0("Saving specific KEGG pathway diagram for: ", target_kegg_id))
pathview(
    gene.data = log2FC_vector,
    pathway.id = target_kegg_id,
    species = "hsa",
    kegg.native = TRUE,
    low = "blue",
    mid = "white",
    high = "red",
    limit = c(min(log2FC_vector), max(log2FC_vector))
)

print("All analyses are complete. Results are saved in the 'clusterprofiler_results' folder.")


# Make sure the directory exists first
# Make sure the directory exists first
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

pathview(
    gene.data = log2FC_vector,
    pathway.id = target_kegg_id,
    species = "hsa",
    kegg.native = TRUE,
    low = "blue",
    mid = "white",
    high = "red",
    limit = c(min(log2FC_vector), max(log2FC_vector)),
    dir = output_dir # ADD THIS LINE to specify the output directory
)
