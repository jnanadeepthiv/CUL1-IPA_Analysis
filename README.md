# INFO B528 - Computational Analysis of High-throughput Biomedical Data
# Project Report — CUL1-IPA Transcriptome Reanalysis

**Author:** Jnana Deepthi Vishnumolakala
**Date:** 04/29/2026

---

## Project Overview

This project presents a transcriptome-wide reanalysis of a publicly available RNA-seq dataset (GSE297241) from Mallick et al. (2026), which characterised CUL1-IPA — a long non-coding RNA derived from intronic polyadenylation at the CUL1 locus.

The original paper demonstrated that CUL1-IPA maintains nucleolar integrity and that its loss disrupts rRNA processing, reduces protein synthesis, and causes G2/M cell cycle arrest, attributing this primarily to MYC suppression. This analysis applies  genome-scale methods to map the full transcriptional consequence of CUL1-IPA loss across three cell lines and identify candidate pathways beyond those the paper examined.

**Key finding:** Interferon alpha signalling is universally activated upon CUL1-IPA loss — a previously unreported downstream consequence confirmed by four independent gene-set databases.

---

## Languages, Tools, and Packages

**Language**
- R (v4.5.2)

**Environment**
- RStudio (v2025.09.2+418)

**Packages**

| Package | Version | Purpose |
|---|---|---|
| limma | 3.66.0 | Differential expression modelling |
| edgeR | 4.8.2 | Expression data utilities |
| ggplot2 | 4.0.2 | Visualisation |
| ggvenn | 0.1.19 | Venn diagrams for DEG overlap |
| pheatmap | 1.0.13 | Heatmaps |
| RColorBrewer | 1.1-3 | Colour palettes |
| ggrepel | 0.9.8 | Label placement in plots |
| fgsea | 1.36.2 | Gene Set Enrichment Analysis (GSEA) |
| msigdbr | 26.1.0 | MSigDB gene set collections (Hallmark, C2) |
| viper | 1.44.0 | Transcription factor activity inference |
| OmnipathR | 3.18.3 | Network resource access for DoRothEA |
| decoupleR | 2.16.0 | TF activity framework |
| clusterProfiler | 4.18.4 | GO and pathway enrichment |
| ReactomePA | 1.54.0 | Reactome pathway analysis |
| org.Hs.eg.db | 3.22.0 | Human gene annotation |
| AnnotationDbi | 1.72.0 | Gene ID conversion |
| dplyr | 1.2.1 | Data manipulation |
| tidyr | 1.3.2 | Data reshaping |
| tibble | 3.3.1 | Data structure utilities |

---

## Workflow - as presented in the course presentation.

1. Load raw expression matrix from GEO (GSE297241) — 18 samples across 3 cell lines
2. Preprocess and filter lowly expressed genes
3. Fit limma-trend differential expression model with cell-line-stratified design matrix, contrasting ASO-treated vs scramble-control within each cell line
4. Extract DEGs (FDR < 0.05) per cell line and identify shared and unique genes using Venn analysis
5. Infer transcription factor activity using VIPER/DoRothEA on the full expression matrix
6. Run GSEA using Hallmark gene sets across all three cell lines
7. Compute ISG activation score and generate ISG heatmap for K562
8. Generate rRNA processing gene heatmap for K562
9. Run C2 immune/stress GSEA for K562

---

## Inputs

**Files:**
- `GSE297241_raw_counts_All_Samples.txt` — expression matrix downloaded from GEO.

**Dataset:**
- 18 samples: 3 cell lines (K562, KMS34, HEK293) × 2 conditions (CUL1-IPA ASO, scramble control) × 3 biological replicates

**Parameters:**
- `FDR < 0.05` — significance threshold for DEGs
- `minSize = 15, maxSize = 500` — GSEA gene set size filters (Hallmark)
- `minSize = 10, maxSize = 500` — GSEA gene set size filters (C2)
- `set.seed(42)` — for GSEA reproducibility

---

## Outputs

All outputs are saved to `Output_limmalog_git/`

**DEG Tables:**
- `DEGs_all/up/down_HEK293/K562/KMS34.csv` — Significant DEGs per cell line
- `Shared_DEGs_all3.csv` — Genes significant in all 3 cell lines
- `Shared_up_all3.csv` — Universally upregulated genes
- `Shared_down_all3.csv` — Universally downregulated genes

**Figures:**
- `Figure1_VIPER_TF_heatmap.png` — TF activity heatmap (top 30 most variable TFs)
- `Figure2_GSEA_Hallmark_HEK293/K562/KMS34.pdf` — GSEA Hallmark barplots
- `Figure3a_ISG_score_K562.pdf` — ISG activation score barplot (K562)
- `Figure3b_ISG_heatmap_K562.png` — Significant ISG heatmap (K562)
- `Figure_rRNA_heatmap_K562.png` — rRNA processing gene heatmap (K562)
- `GSEA_C2_K562.pdf` — C2 immune/stress GSEA barplot (K562)
- `Venn_all/upregulated/downregulated.pdf` — Venn diagrams

**GSEA Tables:**
- `GSEA_Hallmark_HEK293/K562/KMS34.csv` — Full Hallmark results
- `GSEA_C2_K562_significant.csv` — Significant C2 results for K562
- `TF_activity_VIPER_all.csv` — TF activity scores across all cell lines
- `session_info.txt` — R session and package version log

---

## Generative AI Disclosure

Generative AI (ChatGPT) was used as a debugging assistant and for editing assistance in the preparation of the project report.
