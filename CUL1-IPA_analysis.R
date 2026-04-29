# CUL1-IPA RNA-seq Analysis
# Dataset : GSE297241 (18 samples: HEK293, K562, KMS34 x ASO vs scramble)
# Method  : limma

# SECTION 0 — Load libraries

library(OmnipathR)
library(decoupleR)
library(limma)
library(edgeR)
library(ggplot2)
library(ggvenn)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(viper)
library(fgsea)
library(msigdbr)
library(ReactomePA)
library(dplyr)
library(tidyr)
library(tibble)

cat("All libraries loaded.\n")
dir.create("Output_limmalog_git", showWarnings = FALSE)


# SECTION 1 — Load expression matrix

counts_raw <- read.table(
  "input/GSE297241_raw_counts_All_Samples.txt",
  header      = TRUE,
  sep         = "\t",
  check.names = FALSE
)

cat("Raw dimensions:", dim(counts_raw), "\n")

count_cols  <- !colnames(counts_raw) %in% c("gene_id", "symbol")
counts_only <- counts_raw[, count_cols]
rownames(counts_only) <- make.unique(as.character(counts_raw[, "symbol"]))
raw_counts <- as.matrix(counts_only)

cat("Expression matrix dimensions:", dim(raw_counts), "\n")


# SECTION 2 — Sample metadata

# Column order:
#   1-3 : HEK293 ASO    4-6 : HEK293 scramble
#   7-9 : K562 ASO      10-12: K562 scramble
#   13-15: KMS34 ASO   16-18: KMS34 scramble

coldata <- data.frame(
  cell_line = c(rep("HEK293", 6), rep("K562", 6), rep("KMS34", 6)),
  condition = rep(c(rep("ASO", 3), rep("scramble", 3)), 3),
  row.names = colnames(raw_counts)
)

cat("Sample annotation:\n")
print(data.frame(
  sample    = colnames(raw_counts),
  cell_line = coldata$cell_line,
  condition = coldata$condition),
  row.names = FALSE)


# SECTION 3 — Preprocessing and filtering

log2raw_counts <- log2(raw_counts + 1)

cat("\nValue range:", round(range(log2raw_counts), 2), "\n")

keep <- rowSums(log2raw_counts > 1) >= 3
log2raw_counts_filt <- log2raw_counts[keep, ]

cat("Genes before filtering:", nrow(log2raw_counts), "\n")
cat("Genes after filtering: ", nrow(log2raw_counts_filt), "\n")
cat("Genes removed:         ", nrow(log2raw_counts) - nrow(log2raw_counts_filt), "\n\n")


# SECTION 4 — Design matrix and limma model

group       <- factor(paste(coldata$cell_line, coldata$condition, sep = "_"))
design_full <- model.matrix(~ 0 + group)
colnames(design_full) <- levels(group)

cat("Design matrix columns:\n")
print(colnames(design_full))

fit_trend <- lmFit(log2raw_counts_filt, design_full)

contrast_matrix <- makeContrasts(
  HEK293_ASOvsSCR = HEK293_ASO - HEK293_scramble,
  K562_ASOvsSCR   = K562_ASO   - K562_scramble,
  KMS34_ASOvsSCR  = KMS34_ASO  - KMS34_scramble,
  levels = design_full
)

fit_trend2 <- contrasts.fit(fit_trend, contrast_matrix)
fit_trend2 <- eBayes(fit_trend2, trend = TRUE, robust = TRUE)

cat("\nlimma model complete\n")
cat("\nDEG summary (FDR < 0.05):\n")
print(summary(decideTests(fit_trend2, p.value = 0.05)))


# SECTION 5 — Extract DEG results

up_genes   <- function(res) res$gene[!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05 & res$logFC > 0]
down_genes <- function(res) res$gene[!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05 & res$logFC < 0]
sig_genes  <- function(res) res$gene[!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05]

extract_results <- function(fit_obj, coef_name, cell_name) {
  res      <- topTable(fit_obj, coef = coef_name, number = Inf,
                       sort.by = "P", adjust.method = "BH", confint = TRUE)
  res$gene <- rownames(res)
  res$direction <- "NS"
  res$direction[res$adj.P.Val < 0.05 & res$logFC >  0] <- "Up"
  res$direction[res$adj.P.Val < 0.05 & res$logFC <  0] <- "Down"
  write.csv(res, paste0("Output_limmalog_git/limma_full_results_", cell_name, ".csv"), row.names = FALSE)
  sig  <- res[res$adj.P.Val < 0.05, ]
  up   <- sig[sig$logFC > 0, ]
  down <- sig[sig$logFC < 0, ]
  write.csv(sig,  paste0("Output_limmalog_git/DEGs_all_",  cell_name, ".csv"), row.names = FALSE)
  write.csv(up,   paste0("Output_limmalog_git/DEGs_up_",   cell_name, ".csv"), row.names = FALSE)
  write.csv(down, paste0("Output_limmalog_git/DEGs_down_", cell_name, ".csv"), row.names = FALSE)
  cat(cell_name, "Total:", nrow(sig), "Up:", nrow(up), "Down:", nrow(down), "\n")
  return(res)
}

res_HEK293 <- extract_results(fit_trend2, "HEK293_ASOvsSCR", "HEK293")
res_K562   <- extract_results(fit_trend2, "K562_ASOvsSCR",   "K562")
res_KMS34  <- extract_results(fit_trend2, "KMS34_ASOvsSCR",  "KMS34")

v     <- list()
v$E   <- log2raw_counts_filt


# SECTION 6 — Venn diagrams

venn_all  <- list(HEK293 = sig_genes(res_HEK293),
                  K562   = sig_genes(res_K562),
                  KMS34  = sig_genes(res_KMS34))
venn_up   <- list(HEK293 = up_genes(res_HEK293),
                  K562   = up_genes(res_K562),
                  KMS34  = up_genes(res_KMS34))
venn_down <- list(HEK293 = down_genes(res_HEK293),
                  K562   = down_genes(res_K562),
                  KMS34  = down_genes(res_KMS34))

venn_colors <- c("#D85A30", "#534AB7", "#1D9E75")

pdf("Output_limmalog_git/Venn_all.pdf", width = 7, height = 6)
print(ggvenn(venn_all, fill_color = venn_colors, stroke_size = 0.5, set_name_size = 5) +
        labs(title = "All significant DEGs — 3 cell lines (FDR < 0.05)"))
dev.off()

pdf("Output_limmalog_git/Venn_upregulated.pdf", width = 7, height = 6)
print(ggvenn(venn_up, fill_color = venn_colors, stroke_size = 0.5, set_name_size = 5) +
        labs(title = "Upregulated DEGs — gained when CUL1-IPA is lost"))
dev.off()

pdf("Output_limmalog_git/Venn_downregulated.pdf", width = 7, height = 6)
print(ggvenn(venn_down, fill_color = venn_colors, stroke_size = 0.5, set_name_size = 5) +
        labs(title = "Downregulated DEGs — lost when CUL1-IPA is lost"))
dev.off()

shared_all  <- Reduce(intersect, venn_all)
shared_up   <- Reduce(intersect, venn_up)
shared_down <- Reduce(intersect, venn_down)

cat("\nShared all 3:", length(shared_all),
    "| Shared up:", length(shared_up),
    "| Shared down:", length(shared_down), "\n")

write.csv(data.frame(gene = shared_all),  "Output_limmalog_git/Shared_DEGs_all3.csv",  row.names = FALSE)
write.csv(data.frame(gene = shared_up),   "Output_limmalog_git/Shared_up_all3.csv",    row.names = FALSE)
write.csv(data.frame(gene = shared_down), "Output_limmalog_git/Shared_down_all3.csv",  row.names = FALSE)


# SECTION 7 — VIPER TF activity

build_tstat_matrix <- function(res_list, col_names) {
  mat_list <- lapply(res_list, function(res) {
    s <- res$t; s[is.na(s)] <- 0; names(s) <- res$gene; s
  })
  common <- Reduce(intersect, lapply(mat_list, names))
  mat    <- sapply(mat_list, function(s) s[common])
  colnames(mat) <- col_names; rownames(mat) <- common; mat
}

tstat_mat <- build_tstat_matrix(
  list(res_HEK293, res_K562, res_KMS34),
  c("HEK293", "K562", "KMS34")
)

net     <- get_dorothea(organism = "human", levels = c("A", "B", "C"))
tf_acts <- run_viper(mat = tstat_mat, net = net,
                     .source = "source", .target = "target",
                     .mor = "mor", minsize = 4, verbose = FALSE)

tf_matrix <- tf_acts %>%
  select(source, condition, score) %>%
  pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()

write.csv(as.data.frame(tf_matrix),
          "Output_limmalog_git/TF_activity_VIPER_all.csv", row.names = TRUE)

top30 <- names(sort(apply(tf_matrix, 1, var, na.rm = TRUE), decreasing = TRUE))[1:30]

png("Output_limmalog_git/Figure1_VIPER_TF_heatmap.png", width = 1000, height = 1400, res = 150)
pheatmap(tf_matrix[top30, ],
         cluster_cols = FALSE, cluster_rows = TRUE,
         color = colorRampPalette(c("#185FA5", "white", "#993C1D"))(100),
         fontsize_row = 10, fontsize_col = 12,
         main = "TF activity (VIPER/DoRothEA) - top 30 most variable")
dev.off()

cat("\nMYC TF activity:\n"); print(tf_matrix["MYC", ])
cat("\nTop 10 activated TFs:\n"); print(sort(rowMeans(tf_matrix, na.rm = TRUE), decreasing = TRUE)[1:10])
cat("\nTop 10 suppressed TFs:\n"); print(sort(rowMeans(tf_matrix, na.rm = TRUE))[1:10])
cat("Saved: Figure1_VIPER_TF_heatmap.png\n")


# SECTION 8 — GSEA Hallmark

hallmarks <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

make_ranked_list <- function(res) {
  s <- res$t; names(s) <- res$gene
  sort(s[!is.na(s)], decreasing = TRUE)
}

run_gsea_hallmark <- function(res, cell_name) {
  set.seed(42)
  gsea_res <- fgsea(pathways = hallmarks, stats = make_ranked_list(res),
                    minSize = 15, maxSize = 500)
  gsea_res <- gsea_res[order(gsea_res$padj), ]
  write.csv(gsea_res[, c("pathway", "pval", "padj", "NES", "size")],
            paste0("Output_limmalog_git/GSEA_Hallmark_", cell_name, ".csv"), row.names = FALSE)
  
  top20     <- head(gsea_res$pathway, 20)
  plot_data <- gsea_res[gsea_res$pathway %in% top20, ]
  plot_data$short <- gsub("HALLMARK_", "", plot_data$pathway)
  
  p <- ggplot(plot_data, aes(x = NES, y = reorder(short, NES), fill = padj < 0.05)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "#D85A30", "FALSE" = "#B4B2A9"),
                      labels = c("TRUE" = "FDR < 0.05", "FALSE" = "NS")) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(title    = paste("GSEA Hallmark -", cell_name),
         subtitle = "Positive NES = upregulated after CUL1-IPA knockdown",
         x = "NES", y = NULL, fill = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(paste0("Output_limmalog_git/Figure2_GSEA_Hallmark_", cell_name, ".pdf"),
         p, width = 10, height = 8)
  
  sig <- gsea_res[!is.na(gsea_res$padj) & gsea_res$padj < 0.05, ]
  cat(cell_name, "-", nrow(sig), "significant Hallmark pathways\n")
  return(gsea_res)
}

gsea_HEK293 <- run_gsea_hallmark(res_HEK293, "HEK293")
gsea_K562   <- run_gsea_hallmark(res_K562,   "K562")
gsea_KMS34  <- run_gsea_hallmark(res_KMS34,  "KMS34")

cat("\n Key pathway summary \n")
for (gs in c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_G2M_CHECKPOINT",
             "HALLMARK_INTERFERON_ALPHA_RESPONSE",
             "HALLMARK_INTERFERON_GAMMA_RESPONSE")) {
  get_nes <- function(df, pw) {
    r <- df[df$pathway == pw, ]
    if (nrow(r) == 0) return(c("NA", "NA"))
    c(round(r$NES, 2), round(r$padj, 3))
  }
  h <- get_nes(gsea_HEK293, gs)
  k <- get_nes(gsea_K562,   gs)
  m <- get_nes(gsea_KMS34,  gs)
  cat(sprintf("%-40s HEK:%s(p=%s) K562:%s(p=%s) KMS34:%s(p=%s)\n",
              gs, h[1], h[2], k[1], k[2], m[1], m[2]))
}


# SECTION 9 — ISG score and heatmap (K562)

ifn_alpha_genes <- hallmarks[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]]
ifn_gamma_genes <- hallmarks[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]
all_ifn_genes   <- unique(c(ifn_alpha_genes, ifn_gamma_genes))
ifn_in_data     <- intersect(all_ifn_genes, rownames(v$E))

isg_up_HEK293 <- intersect(ifn_in_data, up_genes(res_HEK293))
isg_up_K562   <- intersect(ifn_in_data, up_genes(res_K562))
isg_up_KMS34  <- intersect(ifn_in_data, up_genes(res_KMS34))
isg_shared    <- Reduce(intersect, list(isg_up_HEK293, isg_up_K562, isg_up_KMS34))

cat("\nISGs up HEK293:", length(isg_up_HEK293),
    "K562:", length(isg_up_K562),
    "KMS34:", length(isg_up_KMS34),
    "ALL 3:", length(isg_shared), "\n")
cat("Core ISGs:", paste(isg_shared, collapse = ", "), "\n")

isg_score_genes <- intersect(all_ifn_genes, rownames(v$E))
isg_per_sample  <- colMeans(v$E[isg_score_genes, ])

score_df <- data.frame(
  sample    = names(isg_per_sample),
  isg_score = isg_per_sample,
  cell_line = coldata$cell_line,
  condition = factor(coldata$condition, levels = c("scramble", "ASO"))
)

cat("\nIFN score t-test:\n")
for (cl in c("HEK293", "K562", "KMS34")) {
  sub <- score_df[score_df$cell_line == cl, ]
  tt  <- t.test(isg_score ~ condition, data = sub)
  cat(cl, "p =", round(tt$p.value, 4), "\n")
}

k562_score <- score_df[score_df$cell_line == "K562", ]

p_ifn_k562 <- ggplot(k562_score,
                     aes(x = reorder(sample, isg_score),
                         y = isg_score,
                         fill = condition)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_hline(yintercept = mean(k562_score$isg_score[k562_score$condition == "scramble"]),
             linetype = "dashed", color = "grey40", linewidth = 0.7) +
  scale_fill_manual(values = c("scramble" = "#888780", "ASO" = "#D85A30"),
                    labels = c("scramble" = "Control", "ASO" = "CUL1-IPA KD")) +
  labs(title    = "ISG activation score — K562",
       subtitle = "t-test p = 0.019 | all ASO samples above mean control",
       x = NULL, y = "Mean ISG expression", fill = "Condition") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "top")

ggsave("Output_limmalog_git/Figure3a_ISG_score_K562.pdf", p_ifn_k562, width = 6, height = 5)
cat("Saved: Figure3a_ISG_score_K562.pdf\n")

k562_cols        <- coldata$cell_line == "K562"
mat_isg_k562     <- v$E[ifn_in_data, k562_cols]
mat_isg_k562_sig <- mat_isg_k562[rownames(mat_isg_k562) %in% isg_up_K562, ]

cat("Significant ISGs in K562:", nrow(mat_isg_k562_sig), "\n")

annot_k562_isg <- data.frame(
  Condition = coldata$condition[k562_cols],
  row.names = colnames(mat_isg_k562)
)

core_annotation <- data.frame(
  Core_ISG = ifelse(rownames(mat_isg_k562_sig) %in% isg_shared,
                    "All 3 lines", "K562 only"),
  row.names = rownames(mat_isg_k562_sig)
)

png("Output_limmalog_git/Figure3b_ISG_heatmap_K562.png", width = 800, height = 1000, res = 150)
pheatmap(mat_isg_k562_sig,
         annotation_col    = annot_k562_isg,
         annotation_row    = core_annotation,
         annotation_colors = list(
           Condition = c(ASO = "#E24B4A", scramble = "#888780"),
           Core_ISG  = c("All 3 lines" = "#993C1D", "K562 only" = "#B4B2A9")),
         scale        = "row",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#185FA5", "white", "#993C1D"))(100),
         fontsize_row = 8, fontsize_col = 9,
         main = "Significant ISGs in K562 after CUL1-IPA KD\n(z-scored within K562)")
dev.off()
cat("Saved: Figure3b_ISG_heatmap_K562.png\n")


# SECTION 10 — rRNA heatmap (K562) # highlighted in the paper

rrna_genes <- unique(c("NOP58", "NOP56", "FBL", "DKC1", "GPATCH4",
                       "BOP1", "BRIX1", "GAR1", "NHP2", "NOL9",
                       "NAT10", "NOLC1", "RRP1", "WDR43", "DDX10",
                       "DDX18", "EMG1", "LYAR", "NOB1", "NOL10",
                       "NOL8", "NOP14", "PA2G4", "RRP15", "RRP7A",
                       "TBL3", "TSR1", "UTP14A", "UTP3", "WDR3",
                       "RCL1", "MRTO4", "BYSL"))

rrna_present <- intersect(rrna_genes, rownames(v$E))
cat("rRNA genes found:", length(rrna_present), "\n")

mat_rrna_k562   <- v$E[rrna_present, k562_cols]
annot_k562_rrna <- data.frame(
  Condition = coldata$condition[k562_cols],
  row.names = colnames(mat_rrna_k562)
)

png("Output_limmalog_git/Figure_rRNA_heatmap_K562.png", width = 800, height = 1200, res = 150)
pheatmap(mat_rrna_k562,
         annotation_col    = annot_k562_rrna,
         annotation_colors = list(Condition = c(ASO = "#E24B4A", scramble = "#888780")),
         scale        = "row",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#185FA5", "white", "#993C1D"))(100),
         fontsize_row = 9, fontsize_col = 9,
         main = "rRNA processing genes - K562 (z-scored within cell line)")
dev.off()
cat("Saved: Figure_rRNA_heatmap_K562.png\n")


# SECTION 11 — C2 immune/stress GSEA (K562)

c2_sets <- msigdbr(species = "Homo sapiens", collection = "C2") %>%
  split(x = .$gene_symbol, f = .$gs_name)

search_terms <- paste(
  "INTERFERON|INNATE|STING|CGAS|RIG_I|RIG.I|MDA5|IFIH|",
  "NFKB|NF_KB|NF.KB|TLR|TOLL|AUTOPHAGY|MTOR|RAPAMYCIN|",
  "TP53|P53|MDM2|PKR|EIF2|STRESS_RESPONSE|",
  "NUCLEOL|RIBOSOM|rRNA|APOPTOSIS|ISG|ANTIVIRAL|",
  "CYTOKINE|CHEMOKINE|INTERLEUKIN|JAK|STAT|IFN|TYPE_I|TYPE_II",
  sep = ""
)

immune_stress_names <- grep(search_terms, names(c2_sets), value = TRUE, ignore.case = TRUE)
immune_stress_sets  <- c2_sets[immune_stress_names]
cat("\nC2 immune/stress pathways:", length(immune_stress_names), "\n")

set.seed(42)
gsea_c2_K562 <- fgsea(pathways = immune_stress_sets,
                      stats    = make_ranked_list(res_K562),
                      minSize  = 10, maxSize = 500)
gsea_c2_K562 <- gsea_c2_K562[order(gsea_c2_K562$padj), ]

sig_c2_K562 <- gsea_c2_K562[!is.na(gsea_c2_K562$padj) & gsea_c2_K562$padj < 0.05, ]
write.csv(sig_c2_K562[, c("pathway", "pval", "padj", "NES", "size")],
          "Output_limmalog_git/GSEA_C2_K562_significant.csv", row.names = FALSE)

top20_c2    <- head(gsea_c2_K562$pathway, 20)
plot_c2     <- gsea_c2_K562[gsea_c2_K562$pathway %in% top20_c2, ]
plot_c2$short <- gsub("^KEGG_|^REACTOME_|^BIOCARTA_|^WP_", "", plot_c2$pathway)
plot_c2$short <- substr(plot_c2$short, 1, 55)

p_c2 <- ggplot(plot_c2, aes(x = NES, y = reorder(short, NES), fill = padj < 0.05)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#D85A30", "FALSE" = "#B4B2A9"),
                    labels = c("TRUE" = "FDR < 0.05", "FALSE" = "NS")) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  labs(title = "C2 immune/stress GSEA - K562",
       x = "NES", y = NULL, fill = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 7))

ggsave("Output_limmalog_git/GSEA_C2_K562.pdf", p_c2, width = 12, height = 8)
cat("K562 C2:", nrow(sig_c2_K562), "significant\n")
cat("Saved: GSEA_C2_K562.pdf\n")


# SECTION 12 — Session info

sink("Output_limmalog_git/session_info.txt")
print(sessionInfo())
sink()

pkgs <- c("OmnipathR", "decoupleR", "limma", "edgeR", "ggplot2", "ggvenn",
          "pheatmap", "RColorBrewer", "ggrepel", "clusterProfiler",
          "org.Hs.eg.db", "AnnotationDbi", "viper", "fgsea", "msigdbr",
          "ReactomePA", "dplyr", "tidyr", "tibble")

installed.packages()[pkgs, "Version"]


cat("PIPELINE COMPLETE\n")

cat("Output folder: Output_limmalog_git\n\n")
cat("FIGURES:\n")
cat("  Figure1_VIPER_TF_heatmap.png\n")
cat("  Figure2_GSEA_Hallmark_*.pdf\n")
cat("  Figure3a_ISG_score_K562.pdf\n")
cat("  Figure3b_ISG_heatmap_K562.png\n")
cat("  Figure_rRNA_heatmap_K562.png\n")
cat("  GSEA_C2_K562.pdf\n")
cat("TABLES:\n")
cat("  limma_full_results_*.csv\n")
cat("  DEGs_*.csv\n")
cat("  Shared_DEGs_all3.csv\n")
cat("  GSEA_Hallmark_*.csv\n")
cat("  GSEA_C2_K562_significant.csv\n")
cat("  TF_activity_VIPER_all.csv\n")
cat("  session_info.txt\n")
