# Load necessary RData file
load("~/Allcause_Mortality/allcausemortality.RData") 

# 1-Extract 68 significant metabolic biomarkers in the UK Biobank
library(readxl)
results224 <- read_excel("cox_es_all_224results.xlsx", sheet = 1, na = "NA")

# Load the necessary library
library(dplyr)

# Filter 68 significant metabolic biomarkers with p.adjust < 0.05
significant_metabolites <- results224 %>%
  filter(p.adjust < 0.05) %>%
  select(metabolite)

# Check if there are 68 biomarkers
if(nrow(significant_metabolites) == 68) {
  # Extract the corresponding biomarker variables
  ukb_biomarkers68 <- ukb_all_new_imputed_logall_z %>%
    select(all_of(significant_metabolites$metabolite))
  
  # View the results
  head(ukb_biomarkers68)
} else {
  print(paste("The number of biomarkers after filtering is not 68, actual number is: ", nrow(significant_metabolites)))
}

# 2-WGCNA analysis
# Load necessary packages
library(WGCNA)

# Compute the correlation matrix
correlation_matrix <- cor(ukb_biomarkers68, method = "pearson", use = "pairwise.complete.obs")

# Perform hierarchical clustering to construct a tree
geneTree <- hclust(as.dist(1 - correlation_matrix), method = "average")

# Perform dynamic cutting using cutreeDynamic function
dynamicMods = cutreeDynamic(dendro = geneTree, distM = as.matrix(1 - correlation_matrix),
                            deepSplit = 2, pamRespectsDendro = FALSE)

# Convert module labels to colors
dynamicColors = labels2colors(dynamicMods)

# Print the size of each module
module_sizes <- table(dynamicMods)
print(module_sizes)

# Create a data frame to store biomarkers and their associated module colors
module_biomarkers <- data.frame(
  Biomarker = colnames(ukb_biomarkers68),  # The names of the biomarkers
  Module = dynamicColors  # The color of the module each biomarker belongs to
)

# Save the biomarker-module associations to a CSV file
write.csv(module_biomarkers, file = "ukb_module_biomarkers.csv", row.names = FALSE)

# Sort the correlation matrix and annotation for heatmap
sorted_indices <- order(dynamicColors)
sorted_correlation_matrix <- correlation_matrix[sorted_indices, sorted_indices]
sorted_dynamicColors <- dynamicColors[sorted_indices]
sorted_annotation <- data.frame(Module = factor(sorted_dynamicColors, levels = unique(dynamicColors)))
rownames(sorted_annotation) <- colnames(sorted_correlation_matrix)

# Define custom color range
breaks <- seq(-1, 1, length.out = 100)
color_palette <- colorRampPalette(c("blue", "white", "red"))(99)

# Plot correlation heatmap with white gridlines and transparent background
library(pheatmap)
pdf("ukb_correlation_heatmap.pdf", bg = "transparent", width = 8, height = 7)
pheatmap(sorted_correlation_matrix, 
         clustering_distance_rows = as.dist(1 - sorted_correlation_matrix), 
         clustering_distance_cols = as.dist(1 - sorted_correlation_matrix), 
         clustering_method = "average", 
         annotation_row = sorted_annotation, 
         annotation_col = sorted_annotation, 
         annotation_colors = list(Module = dynamicColors),
         show_rownames = FALSE, 
         show_colnames = FALSE, 
         border_color = "white",  # Add white gridlines
         main = "Correlation Heatmap of Metabolites",
         breaks = breaks,
         color = color_palette)
dev.off()

# Calculate the Topological Overlap Matrix (TOM)
selectTOM <- TOMsimilarityFromExpr(ukb_biomarkers68, power = 6)
dissTOM <- 1 - selectTOM

# Plot the network heatmap
pdf("ukb_network_heatmap.pdf", bg = "transparent", width = 8, height = 7)
TOMplot(dissTOM^7, geneTree, dynamicColors, col = colorRampPalette(c("red", "white", "lightyellow"))(250),
        main = "Network heatmap plot, selected biomarkers")
dev.off()

# Plot the topology map (graphical representation of the TOM)
module_colors <- dynamicColors  # Use the previous module colors

# Plot the heatmap of the TOM matrix
diag(dissTOM) <- NA  # Set diagonal to NA
myheatcol <- colorRampPalette(c("red", "white", "blue"))(250)  # Set color range

# Generate the topology map of the TOM matrix
pdf("ukb_TOM_topology.pdf", bg = "transparent", width = 8, height = 7)
plotTOM = pheatmap(dissTOM, 
                   color = myheatcol,
                   show_rownames = FALSE, 
                   show_colnames = FALSE, 
                   cluster_rows = geneTree, 
                   cluster_cols = geneTree,
                   annotation_row = data.frame(Module = module_colors),
                   annotation_colors = list(Module = module_colors),
                   main = "TOM Topology Heatmap")
dev.off()

# Save the TOM matrix and dynamic cut results
write.csv(dissTOM, file = "ukb_dissTOM_matrix.csv", row.names = TRUE)
write.csv(module_sizes, file = "ukb_module_sizes.csv", row.names = TRUE)
