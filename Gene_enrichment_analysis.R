###################This script was used for performing gene enrichment analysis for the bagworm genes with evolutionary signatures (using Fisher's exact test)

library(data.table)
kegg <- fread("KEGG_pathways_annotation_list.tsv", header = FALSE)
colnames(kegg) <- c("gene", "ko")
focal <- scan("focal_gene_list.txt", what = "", quiet = TRUE)
background <- scan("background_gene_list.txt", what = "", quiet = TRUE)
all_genes <- union(focal, background)
kegg_all <- kegg[gene %in% all_genes]
kegg_all[, is_focal := ifelse(gene %in% focal, 1, 0)]
counts <- kegg_all[, .(focal = sum(is_focal), background = .N), by = ko]
total_focal <- length(unique(focal))
total_background <- length(unique(all_genes))
counts[, fold_enrichment := (focal / total_focal) / (background / total_background)]

counts[, pvalue := apply(.SD, 1, function(x) {
  a <- x["focal"]
  b <- total_focal - a
  c <- x["background"] - a  # background-only count
  d <- total_background - x["background"]
  fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
}), .SDcols = c("focal", "background")]

counts[, padj := p.adjust(pvalue, method = "BH")]
significant <- counts[padj < 0.05]
significant <- significant[order(-focal)]
fwrite(significant, "KEGG_enrichment_result.tsv", sep = "\t")

###################Plotting the results

library(data.table)
library(ggplot2)
df <- fread("KEGG_enrichment_result.tsv")
df <- df[padj < 0.01]
df <- df[order(padj)][1:20]
df[, neglogpadj := -log10(padj)]
df[, ko := factor(ko, levels = rev(ko))]

pdf(file = 'Gene_enrichment.pdf', width=10, height=10)
ggplot(df, aes(x = fold_enrichment, y = ko)) +
  geom_point(aes(size = focal, color = neglogpadj)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "KEGG Pathway Enrichment",
       x = "Fold Enrichment",
       y = "KEGG Pathway",
       color = "-log10(FDR)",
       size = "Focal gene count") +
  theme_gray(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  )
dev.off()
