# Proteomics Differential Expression + Enrichment Pipeline (Mouse) — 6 Contrasts

This repository contains an **R-based analysis pipeline** for DIA proteomics results exported from Excel. It generates:

- **Interactive volcano plots** (HTML + PNG) per contrast
- **Differentially expressed proteins (DEPs/DEGs)** tables per contrast (IPA/STRING-ready)
- **Top variable DEPs heatmap** across contrasts (logFC-based)
- **Over-representation analysis (ORA)** for **GO Biological Process** and **Reactome**
- **Gene Set Enrichment Analysis (GSEA)** for **GO Biological Process** and **Reactome**, including a **single combined plot** across all contrasts
- **DEG overlap summary** (UpSet plot + counts table)

> Organism: **Mus musculus (mouse)**  
> Annotation mapping: **UniProt → Entrez** via `org.Mm.eg.db`

---
01_volcano/ # Volcano plots (HTML + PNG) per contrast
02_DEG_tables_IPA/ # FULL results + DEG-only tables + IPA/STRING lists
03_heatmap/ # Heatmap image
04_ORA_enrichment/ # ORA GO + Reactome results and combined plots
05_GSEA_enrichment/ # GSEA GO + Reactome results and combined plots
06_DEG_overlap/ # UpSet overlap plot + DEG counts table
## Project Structure

After running the pipeline, the following folders are created automatically:

---

## Input Data

### Required file
Place the Excel file in the project root:

- `CompadreCM_20250520_01_DIA_results.xlsx`

### Required sheet
The pipeline reads:

- Sheet: `Protein Results`
- Skips the first `2` rows (`skip = 2`)

### Required columns
Your Excel file must include:

- A protein header column (one of):
  - `Protein Annotation` **(preferred)**
  - `Protein_Annotation`
  - otherwise the **first column** is used as a fallback

- LogFC columns for each contrast (exact names must match):
  - `Tocoflexol_IR_vs_Vehicle_IR`
  - `Tocoflexol_vs_Vehicle`
  - `Vehicle_IR_vs_Vehicle`
  - `Tocoflexol_IR_vs_Tocoflexol`
  - `Tocoflexol_IR_vs_Vehicle`
  - `Tocoflexol_vs_Vehicle_IR`

- Adjusted p-value columns (padj):  
  The pipeline assumes **padj is located 6 columns to the right** of each logFC column (based on the Excel export format used).

If your file uses a different layout for adjusted p-values, update the `offset = 6` setting in the script function `get_padj_col()`.

---

## Installation / Requirements

### R version
- R (>= 4.1 recommended)

### Packages
The script will auto-install missing packages using CRAN and Bioconductor.

CRAN packages include:
- `readxl`, `dplyr`, `ggplot2`, `plotly`, `htmlwidgets`, `webshot`, `pheatmap`, `tidyr`, `stringr`, `UpSetR`

Bioconductor packages include:
- `clusterProfiler`, `org.Mm.eg.db`, `ReactomePA`, `AnnotationDbi`

### PhantomJS (for saving Plotly → PNG)
This pipeline uses `webshot`, which relies on **PhantomJS**.  
The script runs:

```r
webshot::install_phantomjs()

## How to Run
Put your Excel file in the repo root:

CompadreCM_20250520_01_DIA_results.xlsx

Open RStudio and run:

source("analysis_pipeline.R")

Outputs will be written into the folders listed in Project Structure.

Key Parameters You Can Edit

At the top of the script under 0) CONFIG:

DEG thresholds (used for volcano + DEG tables + ORA)

deg_padj_cutoff (default: 0.05)

deg_logfc_cutoff (default: 1)

ORA thresholds

ora_pvalue_cutoff (default: 0.05)

ora_min_mapped (default: 5)

GSEA thresholds (more permissive by design)

gsea_pvalue_cutoff (default: 0.2)

gsea_prefilter_logfc (default: 0.2)

gsea_prefilter_padj (default: 0.2)

Heatmap

topN_heatmap (default: 50)

Outputs Explained
Volcano Plots (01_volcano/)

Interactive .html and static .png per contrast

Points colored by significance:

Upregulated (red)

Downregulated (blue)

Not significant (grey)

DEG Tables (02_DEG_tables_IPA/)

For each contrast:

*_FULL_results.csv (all proteins)

*_DEG_Table_IPA.csv (DEGs only) with columns:

UniProtID (e.g., A0A075B5J9)

UniProtName (e.g., A0A075B5J9_MOUSE)

ProteinName (descriptive protein name)

logFC, padj, Regulation

*_IPA_STRING_list.txt (simple list for IPA/STRING)

Heatmap (03_heatmap/)

Top_DEPs_Heatmap_labeled.png

Built from DEG-only proteins across contrasts using logFC

Rows are labeled as:

UniProtID | ProteinName

ORA Enrichment (04_ORA_enrichment/)

Per contrast:

*_ORA_GO.csv

*_ORA_Reactome.csv

Combined:

ORA_GO_All_Contrasts.csv

ORA_Reactome_All_Contrasts.csv

Combined barplots for GO and Reactome

GSEA Enrichment (05_GSEA_enrichment/)

Combined tables:

GSEA_GO_All_Contrasts.csv

GSEA_Reactome_All_Contrasts.csv

Plots:

Facet plots (top 10 terms per contrast)

Single combined plots (global top terms across all contrasts)

Interpretation:

NES > 0 = enriched in the “up” direction

NES < 0 = enriched in the “down” direction

DEG Overlap (06_DEG_overlap/)

DEG_Overlap_UpSet.pdf

DEG_summary_counts.csv

Notes / Troubleshooting
1) “padj column not found”

Your Excel format might not place padj 6 columns after logFC.
Fix by updating the offset in get_padj_col().

2) Mapping failures (UniProt → Entrez)

Some UniProt IDs may not map (common with fragments/unreviewed entries).
The pipeline removes unmapped IDs automatically.

3) Empty enrichment results

If a contrast has few DEGs or weak signal, ORA/GSEA may return no significant terms under the chosen cutoff.
Try:

relaxing DEG thresholds (deg_padj_cutoff, deg_logfc_cutoff)

relaxing GSEA cutoffs (gsea_pvalue_cutoff)

## Citation / Acknowledgements

This analysis uses:

clusterProfiler for GO and GSEA

ReactomePA for Reactome ORA/GSEA

org.Mm.eg.db for mouse gene annotations

## Author

Pipeline assembled and customized for DIA proteomics multi-contrast analysis.
