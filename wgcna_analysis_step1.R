#!/usr/bin/env Rscript
#
# WGCNA Step 1: Module Detection + Module-Trait Correlation + Hub Genes
# Adapts Jia et al. (2026) J Anim Sci Biotechnol framework for pig study
#
# Usage: Rscript wgcna_analysis_step1.R <tissue> <expr_file> <trait_file> <out_prefix>
#   tissue: "Liver" or "Muscle"
#

args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
expr_file <- args[2]
trait_file <- args[3]
out_prefix <- args[4]

suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads(4)

cat(sprintf("\n========================================\n"))
cat(sprintf("WGCNA Analysis: %s\n", tissue))
cat(sprintf("========================================\n"))

# ---- 1. Load data ----
cat("Loading expression data...\n")
expr_data <- read.csv(expr_file, row.names = 1, check.names = FALSE)
cat(sprintf("  Samples: %d, Genes: %d\n", nrow(expr_data), ncol(expr_data)))

trait_data <- read.csv(trait_file, row.names = 1, check.names = FALSE)
cat(sprintf("  Traits: %s\n", paste(colnames(trait_data), collapse = ", ")))

# CSV is already (samples × genes) — WGCNA expects samples as rows, genes as columns
datExpr <- expr_data

# Check for missing values and genes with zero variance
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    cat(sprintf("  Removing %d problematic genes\n", sum(!gsg$goodGenes)))
    datExpr <- datExpr[, gsg$goodGenes]
  }
}

# ---- 2. Soft-thresholding power selection ----
cat("\nSelecting soft-thresholding power...\n")
powers <- c(1:20, seq(22, 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 2, networkType = "signed")

pdf(sprintf("%s_soft_threshold.pdf", out_prefix), width = 8, height = 4)
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = paste("Scale independence (", tissue, ")", sep = ""))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.8, col = "blue", lty = 2)

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = paste("Mean connectivity (", tissue, ")", sep = ""))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# Select power: first power where R^2 > 0.8, with minimum fallback
soft_power <- sft$powerEstimate
if (is.na(soft_power)) {
  # Find closest to R^2 = 0.8
  idx <- which.min(abs(-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2] - 0.8))
  soft_power <- sft$fitIndices[idx, 1]
}
# Enforce minimum soft power for signed networks (too low = no co-expression structure)
if (soft_power < 6) {
  cat(sprintf("  Power %d too low, forcing soft_power = 6 for adequate co-expression structure\n", soft_power))
  soft_power <- 6
}
cat(sprintf("  Selected soft power: %d\n", soft_power))

# ---- 3. Network construction and module detection ----
cat("\nConstructing co-expression network...\n")

# For large datasets, use blockwiseModules
n_genes <- ncol(datExpr)

if (n_genes > 20000) {
  cat("  Using blockwiseModules (large dataset)\n")
  # First filter to most variable genes for practical runtime
  var_expr <- apply(datExpr, 2, var)
  keep_idx <- var_expr >= quantile(var_expr, 0.25)  # Keep top 75% most variable
  datExpr_filtered <- datExpr[, keep_idx]
  cat(sprintf("  Filtered to %d most variable genes\n", ncol(datExpr_filtered)))

  net <- blockwiseModules(datExpr_filtered, power = soft_power,
                          TOMType = "signed", minModuleSize = 20,
                          reassignThreshold = 0, mergeCutHeight = 0.15,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = FALSE, verbose = 2, maxBlockSize = 15000,
	                          deepSplit = 3)
  # Map back to full gene set
  gene_filtered <- colnames(datExpr_filtered)
} else {
  cat("  Using standard one-step network construction\n")
  adjacency <- adjacency(datExpr, power = soft_power, type = "signed")
  TOM <- TOMsimilarity(adjacency)
  dissTOM <- 1 - TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average")

  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit = 3, pamRespectsDendro = FALSE,
                               minClusterSize = 20)
  dynamicColors <- labels2colors(dynamicMods)

  MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
  MEs <- MEList$eigengenes
  MEDiss <- 1 - cor(MEs)
  METree <- hclust(as.dist(MEDiss), method = "average")

  merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.15, verbose = 3)
  mergedColors <- merge$colors
  mergedMEs <- merge$newMEs

  net <- list(colors = mergedColors, MEs = mergedMEs, dendrograms = geneTree)
  gene_filtered <- colnames(datExpr)
}

cat(sprintf("  Modules detected: %d\n", length(unique(net$colors))))

# ---- 4. Module-Trait Correlation ----
cat("\nComputing module-trait correlations...\n")

# Get module eigengenes
MEs <- net$MEs
if (is.null(MEs)) {
  MEs <- moduleEigengenes(datExpr_filtered, colors = net$colors)$eigengenes
}

# Remove grey module (unassigned genes)
module_colors <- net$colors
if (is.list(module_colors)) module_colors <- module_colors[[1]]
if (is.null(names(module_colors))) names(module_colors) <- gene_filtered

unique_modules <- setdiff(unique(module_colors), "grey")
ME_cols <- paste("ME", unique_modules, sep = "")
MEs <- MEs[, ME_cols, drop = FALSE]
colnames(MEs) <- unique_modules

# Align trait data with expression samples
common_samples <- intersect(rownames(datExpr), rownames(trait_data))
cat(sprintf("  Common samples: %d\n", length(common_samples)))

MEs_aligned <- MEs[common_samples, , drop = FALSE]
trait_aligned <- trait_data[common_samples, , drop = FALSE]

# Compute correlations
module_trait_cor <- cor(MEs_aligned, trait_aligned, use = "pairwise.complete.obs")
module_trait_pvalue <- corPvalueStudent(module_trait_cor, nrow(MEs_aligned))

# ---- 5. Save module-trait heatmap ----
pdf(sprintf("%s_module_trait_heatmap.pdf", out_prefix), width = 10, height = max(8, ncol(MEs) * 0.4))
labeledHeatmap(Matrix = module_trait_cor,
               xLabels = colnames(trait_aligned),
               yLabels = rownames(module_trait_cor),
               ySymbols = rownames(module_trait_cor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = round(module_trait_cor, 2),
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               main = paste(tissue, "Module-Trait Correlations"))
dev.off()

# ---- 6. Gene Significance and Module Membership ----
cat("\nComputing gene significance and module membership...\n")

# For each trait, compute GS (correlation of gene with trait)
# and for each module, compute MM (kME)
gs_all <- list()
for (trait_name in colnames(trait_aligned)) {
  trait_vec <- trait_aligned[common_samples, trait_name]
  gs <- as.numeric(cor(datExpr[common_samples, ], trait_vec, use = "pairwise.complete.obs"))
  names(gs) <- colnames(datExpr)
  gs_all[[trait_name]] <- gs
}

# Module membership for all genes
kME_all <- signedKME(datExpr[common_samples, ], MEs_aligned)
colnames(kME_all) <- paste("kME", unique_modules, sep = "_")

# Compile gene-level results
gene_info <- data.frame(
  Gene = colnames(datExpr),
  Module = ifelse(colnames(datExpr) %in% gene_filtered,
                  module_colors[match(colnames(datExpr), gene_filtered)], "grey"),
  stringsAsFactors = FALSE
)

for (tn in names(gs_all)) {
  gene_info[[paste("GS", tn, sep = "_")]] <- gs_all[[tn]][gene_info$Gene]
}

# Add kME for the gene's own module
gene_info$kME_module <- NA
for (i in 1:nrow(gene_info)) {
  mod <- gene_info$Module[i]
  if (mod != "grey") {
    kme_col <- paste("kME", mod, sep = "_")
    if (kme_col %in% colnames(kME_all)) {
      gene_info$kME_module[i] <- kME_all[gene_info$Gene[i], kme_col]
    }
  }
}

# ---- 7. Hub Gene Identification ----
cat("\nIdentifying hub genes...\n")

hub_genes <- data.frame()
for (mod in unique_modules) {
  mod_genes <- gene_info[gene_info$Module == mod, ]
  if (nrow(mod_genes) < 10) next

  # Sort by kME (module membership) within module
  mod_genes <- mod_genes[order(-mod_genes$kME_module), ]

  # Top 20 hub genes per module
  n_hub <- min(20, nrow(mod_genes))
  top_hubs <- mod_genes[1:n_hub, ]

  # Mean trait correlation for module
  mod_cor_row <- module_trait_cor[mod, , drop = FALSE]

  hub_genes <- rbind(hub_genes, data.frame(
    Module = mod,
    Module_Size = nrow(mod_genes),
    Gene = top_hubs$Gene,
    kME = top_hubs$kME_module,
    top_hubs[, grep("^GS_", colnames(top_hubs)), drop = FALSE],
    stringsAsFactors = FALSE
  ))
}

# ---- 8. Save Results ----
cat("\nSaving results...\n")

# Module assignment
write.csv(gene_info, sprintf("%s_gene_module_assignment.csv", out_prefix), row.names = FALSE)

# Module-trait correlations
write.csv(module_trait_cor, sprintf("%s_module_trait_cor.csv", out_prefix))
write.csv(module_trait_pvalue, sprintf("%s_module_trait_pvalue.csv", out_prefix))

# Hub genes
write.csv(hub_genes, sprintf("%s_hub_genes.csv", out_prefix), row.names = FALSE)

# Module sizes
module_sizes <- table(gene_info$Module)
write.csv(data.frame(Module = names(module_sizes), Size = as.numeric(module_sizes)),
          sprintf("%s_module_sizes.csv", out_prefix), row.names = FALSE)

# Save key parameters
cat(sprintf("SoftPower: %d\n", soft_power))
cat(sprintf("Total modules: %d\n", length(unique_modules)))
cat(sprintf("Genes assigned to modules: %d / %d\n",
            sum(gene_info$Module != "grey"), nrow(gene_info)))

cat(sprintf("\n%s WGCNA complete!\n", tissue))
cat(sprintf("Output files: %s_*\n", out_prefix))
