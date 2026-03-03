############################################################
# Proteomics 6-contrast analysis pipeline
# - Interactive volcano plots (HTML + PNG)
# - DEG tables (IPA/STRING-ready) with UniProtID, UniProtName, ProteinName
# - Heatmap (Top-N variable DEPs across contrasts, using DEG-only proteins)
# - ORA enrichment: GO BP + Reactome (barplots + CSV exports)
# - GSEA enrichment: GO BP + Reactome (per-contrast + combined “single plot” + CSV exports)
# - DEG overlap summary (UpSet + counts)
############################################################

############################
# 0) CONFIG
############################
# Input
file_path <- "CompadreCM_20250520_01_DIA_results.xlsx"
sheet_name <- "Protein Results"
skip_rows <- 2

# Contrasts (must match columns present in the Excel file)
contrasts <- c(
  "Tocoflexol_IR_vs_Vehicle_IR",
  "Tocoflexol_vs_Vehicle",
  "Vehicle_IR_vs_Vehicle",
  "Tocoflexol_IR_vs_Tocoflexol",
  "Tocoflexol_IR_vs_Vehicle",
  "Tocoflexol_vs_Vehicle_IR"
)

# Volcano/DEG thresholds (adjust as you like)
deg_padj_cutoff  <- 0.05
deg_logfc_cutoff <- 1

# ORA thresholds
ora_pvalue_cutoff <- 0.05
ora_min_mapped    <- 5

# GSEA thresholds
gsea_pvalue_cutoff <- 0.2
gsea_min_gs_size   <- 5
gsea_max_gs_size   <- 500
# Optional prefilter for ranked list construction (keeps more genes than DEGs)
gsea_prefilter_logfc <- 0.2
gsea_prefilter_padj  <- 0.2
gsea_min_ranked_genes <- 10

# Heatmap settings
topN_heatmap <- 50

# Outputs
out_dir_volcano   <- "01_volcano"
out_dir_deg       <- "02_DEG_tables_IPA"
out_dir_heatmap   <- "03_heatmap"
out_dir_ora       <- "04_ORA_enrichment"
out_dir_gsea      <- "05_GSEA_enrichment"
out_dir_overlap   <- "06_DEG_overlap"

dir.create(out_dir_volcano, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_deg,     showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_heatmap, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_ora,     showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_gsea,    showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_overlap, showWarnings = FALSE, recursive = TRUE)


############################
# 1) PACKAGES
############################
pkgs_cran <- c("readxl", "dplyr", "ggplot2", "plotly", "htmlwidgets", "webshot",
               "pheatmap", "tidyr", "stringr")
pkgs_bioc <- c("clusterProfiler", "org.Mm.eg.db", "ReactomePA", "AnnotationDbi")
pkgs_extra <- c("UpSetR")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

for (p in pkgs_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in pkgs_extra) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in pkgs_bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)

library(readxl)
library(dplyr)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(webshot)
library(pheatmap)
library(tidyr)
library(stringr)

library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(AnnotationDbi)

library(UpSetR)

# PhantomJS for webshot (run once; safe to keep)
webshot::install_phantomjs()


############################
# 2) LOAD DATA
############################
df <- readxl::read_excel(file_path, sheet = sheet_name, skip = skip_rows)

# Figure out which protein annotation column exists
protein_col <- dplyr::case_when(
  "Protein Annotation"  %in% colnames(df) ~ "Protein Annotation",
  "Protein_Annotation"  %in% colnames(df) ~ "Protein_Annotation",
  TRUE ~ colnames(df)[1]  # fallback: first column
)

# Helper: find padj column for a given contrast
# Your file seems to store padj 6 columns after logFC (based on earlier code).
get_padj_col <- function(contrast_name, df_cols, offset = 6) {
  idx <- which(df_cols == contrast_name)
  if (length(idx) == 0) return(NA_character_)
  j <- idx + offset
  if (j > length(df_cols)) return(NA_character_)
  df_cols[j]
}

# Helper: parse FASTA-like headers into UniProt fields
parse_protein_fields <- function(protein_header) {
  # Examples:
  # >tr|A0A075B5J9|A0A075B5J9_MOUSE Immunoglobulin ... OS=Mus musculus ...
  uniprot_id   <- sub("^>\\w+\\|(\\w+)\\|.*$", "\\1", protein_header)
  uniprot_name <- sub("^>\\w+\\|\\w+\\|(\\S+).*$", "\\1", protein_header)     # A0A..._MOUSE
  prot_name    <- sub("^>\\w+\\|\\w+\\|\\S+\\s+", "", protein_header)         # drop header part
  prot_name    <- sub(" OS=.*$", "", prot_name)                               # drop OS=...
  tibble(
    UniProtID = uniprot_id,
    UniProtName = uniprot_name,
    ProteinName = prot_name
  )
}


############################
# 3) VOLCANO PLOT (INTERACTIVE) FUNCTION
############################
plot_volcano_plotly <- function(data, contrast_name, out_dir,
                                pcut = 0.05, logfccut = 1) {
  data <- data %>%
    mutate(
      Sig = case_when(
        padj < pcut & logFC >  logfccut ~ "Upregulated",
        padj < pcut & logFC < -logfccut ~ "Downregulated",
        TRUE                             ~ "Not significant"
      ),
      negLogP = -log10(padj)
    )
  
  color_map <- c(
    "Downregulated"   = "blue",
    "Upregulated"     = "red",
    "Not significant" = "grey"
  )
  
  p <- ggplot(data, aes(
    x = logFC,
    y = negLogP,
    color = Sig,
    text = paste0(
      "Protein: ", Protein, "<br>",
      "logFC: ", round(logFC, 3), "<br>",
      "adj.P: ", signif(padj, 3), "<br>",
      "Class: ", Sig
    )
  )) +
    geom_point(alpha = 0.7, size = 2.2) +
    scale_color_manual(values = color_map) +
    geom_vline(xintercept = c(-logfccut, logfccut), linetype = "dashed") +
    geom_hline(yintercept = -log10(pcut), linetype = "dashed") +
    theme_minimal() +
    labs(
      title = paste("Interactive Volcano Plot:", contrast_name),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Regulation"
    )
  
  p_interactive <- ggplotly(p, tooltip = "text")
  
  html_file <- file.path(out_dir, paste0(contrast_name, "_volcano_plotly.html"))
  png_file  <- file.path(out_dir, paste0(contrast_name, "_volcano_plotly.png"))
  
  htmlwidgets::saveWidget(p_interactive, html_file, selfcontained = TRUE)
  webshot::webshot(html_file, file = png_file, vwidth = 1200, vheight = 900, delay = 1)
  
  invisible(p_interactive)
}


############################
# 4) BUILD PER-CONTRAST TABLES + DEG TABLE EXPORT (IPA/STRING-READY)
############################
deg_list <- list()          # stores DEG protein headers (full string) per contrast
deg_uniprot_list <- list()  # stores DEG UniProtIDs per contrast (for ORA mapping)
contrast_tables <- list()   # stores df_contrast per contrast

for (contrast in contrasts) {
  logfc_col <- contrast
  padj_col  <- get_padj_col(contrast, colnames(df), offset = 6)
  
  if (!(logfc_col %in% colnames(df))) {
    message("Skipping (missing logFC col): ", contrast)
    next
  }
  if (is.na(padj_col) || !(padj_col %in% colnames(df))) {
    message("Skipping (missing padj col): ", contrast)
    next
  }
  
  df_contrast <- df %>%
    transmute(
      Protein = .data[[protein_col]],
      logFC   = suppressWarnings(as.numeric(.data[[logfc_col]])),
      padj    = suppressWarnings(as.numeric(.data[[padj_col]]))
    ) %>%
    filter(!is.na(logFC) & !is.na(padj))
  
  # Add UniProtID / UniProtName / ProteinName
  parsed <- parse_protein_fields(df_contrast$Protein)
  df_contrast <- bind_cols(df_contrast, parsed)
  
  # DEG classification
  df_contrast <- df_contrast %>%
    mutate(
      Regulation = case_when(
        padj < deg_padj_cutoff & logFC >  deg_logfc_cutoff ~ "Upregulated",
        padj < deg_padj_cutoff & logFC < -deg_logfc_cutoff ~ "Downregulated",
        TRUE                                               ~ "Not significant"
      )
    )
  
  # Save full results table for this contrast
  write.csv(df_contrast,
            file = file.path(out_dir_deg, paste0(contrast, "_FULL_results.csv")),
            row.names = FALSE)
  
  # Save DEG-only table (IPA-friendly)
  deg_df <- df_contrast %>% filter(Regulation != "Not significant")
  
  write.csv(
    deg_df %>% select(UniProtID, UniProtName, ProteinName, logFC, padj, Regulation),
    file = file.path(out_dir_deg, paste0(contrast, "_DEG_Table_IPA.csv")),
    row.names = FALSE
  )
  
  # IPA/STRING simple list (UniProtName is often convenient; switch to UniProtID if you prefer)
  write.table(
    deg_df$UniProtName,
    file = file.path(out_dir_deg, paste0(contrast, "_IPA_STRING_list.txt")),
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  
  # Save to lists for downstream steps
  deg_list[[contrast]]        <- deg_df$Protein
  deg_uniprot_list[[contrast]] <- deg_df$UniProtID
  contrast_tables[[contrast]] <- df_contrast
  
  # Volcano plot for this contrast
  plot_volcano_plotly(df_contrast %>% select(Protein, logFC, padj),
                      contrast_name = contrast,
                      out_dir = out_dir_volcano,
                      pcut = deg_padj_cutoff,
                      logfccut = deg_logfc_cutoff)
  
  message("Done: ", contrast, " | DEGs: ", nrow(deg_df))
}

cat("✅ Volcano plots + DEG tables exported.\n")


############################
# 5) HEATMAP: Top-N variable DEPs across contrasts (DEG-only proteins)
############################
# Build logFC matrix using DEG-only proteins per contrast
logfc_matrix_list <- list()

for (contrast in names(deg_list)) {
  logfc_col <- contrast
  if (!(logfc_col %in% colnames(df))) next
  
  degs <- deg_list[[contrast]]
  tmp <- df %>%
    transmute(
      Protein = .data[[protein_col]],
      logFC   = suppressWarnings(as.numeric(.data[[logfc_col]]))
    ) %>%
    filter(Protein %in% degs)
  
  colnames(tmp) <- c("Protein", contrast)
  logfc_matrix_list[[contrast]] <- tmp
}

if (length(logfc_matrix_list) >= 2) {
  logfc_combined <- Reduce(function(x, y) merge(x, y, by = "Protein", all = TRUE), logfc_matrix_list)
  rownames(logfc_combined) <- logfc_combined$Protein
  logfc_matrix <- as.matrix(logfc_combined[, setdiff(colnames(logfc_combined), "Protein"), drop = FALSE])
  
  # row-scale
  logfc_matrix_scaled <- t(scale(t(logfc_matrix)))
  logfc_matrix_scaled[!is.finite(logfc_matrix_scaled)] <- 0
  
  # top-N by variance (guard against out-of-bounds)
  row_vars <- apply(logfc_matrix_scaled, 1, var, na.rm = TRUE)
  row_vars <- row_vars[is.finite(row_vars)]
  topN <- min(topN_heatmap, length(row_vars))
  top_proteins <- names(sort(row_vars, decreasing = TRUE))[seq_len(topN)]
  top_matrix <- logfc_matrix_scaled[top_proteins, , drop = FALSE]
  
  # Label rows as "UniProtID | ProteinName"
  row_labels <- sapply(rownames(top_matrix), function(h) {
    id   <- sub("^>\\w+\\|(\\w+)\\|.*$", "\\1", h)
    name <- sub("^>\\w+\\|\\w+\\|\\S+\\s+", "", h)
    name <- sub(" OS=.*$", "", name)
    paste(id, "|", name)
  })
  rownames(top_matrix) <- row_labels
  
  # Column annotation (edit as needed)
  # Must have rownames exactly matching the contrast column names used in the matrix
  contrast_annotation <- data.frame(
    Treatment  = c("Tocoflexol", "Tocoflexol", "Vehicle", "Tocoflexol", "Tocoflexol", "Tocoflexol"),
    Irradiated = c("Yes",        "No",         "Yes",     "Yes",        "Yes",        "No")
  )
  rownames(contrast_annotation) <- contrasts
  
  # Keep only annotation rows for columns that exist in matrix
  contrast_annotation <- contrast_annotation[intersect(colnames(top_matrix), rownames(contrast_annotation)), , drop = FALSE]
  
  png(file.path(out_dir_heatmap, "Top_DEPs_Heatmap_labeled.png"),
      width = 3000, height = 2000, res = 300)
  pheatmap(
    top_matrix,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = contrast_annotation,
    color = colorRampPalette(c("blue", "yellow", "red"))(100),
    main = paste0("Top ", topN, " Most Variable DEPs (logFC) across contrasts"),
    fontsize_row = 6,
    fontsize_col = 11,
    border_color = NA
  )
  dev.off()
  
  cat("✅ Heatmap saved.\n")
} else {
  message("⚠️ Not enough contrasts with DEG logFC values to build heatmap.")
}


############################
# 6) ORA ENRICHMENT (DEG LISTS): GO BP + Reactome
############################
ora_go_list <- list()
ora_reactome_list <- list()

for (contrast in names(deg_uniprot_list)) {
  message("ORA enrichment: ", contrast)
  
  unip <- deg_uniprot_list[[contrast]]
  unip <- unique(na.omit(unip))
  
  if (length(unip) < ora_min_mapped) {
    message("  ⚠️ Too few UniProt IDs, skipping ORA: ", length(unip))
    next
  }
  
  entrez_ids <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = unip,
    column = "ENTREZID",
    keytype = "UNIPROT",
    multiVals = "first"
  )
  entrez_ids <- unique(na.omit(entrez_ids))
  
  if (length(entrez_ids) < ora_min_mapped) {
    message("  ⚠️ Too few mapped Entrez IDs, skipping ORA: ", length(entrez_ids))
    next
  }
  
  go_res <- clusterProfiler::enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pvalueCutoff  = ora_pvalue_cutoff,
    readable      = TRUE
  )
  
  re_res <- ReactomePA::enrichPathway(
    gene         = entrez_ids,
    organism     = "mouse",
    pvalueCutoff = ora_pvalue_cutoff,
    readable     = TRUE
  )
  
  go_df <- as.data.frame(go_res)
  re_df <- as.data.frame(re_res)
  
  if (nrow(go_df) > 0) {
    go_df$Contrast <- contrast
    ora_go_list[[contrast]] <- go_df
    write.csv(go_df, file.path(out_dir_ora, paste0(contrast, "_ORA_GO.csv")), row.names = FALSE)
  }
  if (nrow(re_df) > 0) {
    re_df$Contrast <- contrast
    ora_reactome_list[[contrast]] <- re_df
    write.csv(re_df, file.path(out_dir_ora, paste0(contrast, "_ORA_Reactome.csv")), row.names = FALSE)
  }
}

ora_go <- bind_rows(ora_go_list)
ora_reactome <- bind_rows(ora_reactome_list)

write.csv(ora_go,       file.path(out_dir_ora, "ORA_GO_All_Contrasts.csv"), row.names = FALSE)
write.csv(ora_reactome, file.path(out_dir_ora, "ORA_Reactome_All_Contrasts.csv"), row.names = FALSE)

# Combined ORA barplots (top 10 per contrast by p.adjust)
plot_ora_bar <- function(df, title, outfile, top_n = 10) {
  if (nrow(df) == 0) return(invisible(NULL))
  
  top_df <- df %>%
    group_by(Contrast) %>%
    slice_min(order_by = p.adjust, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(EnrichmentScore = -log10(p.adjust))
  
  p <- ggplot(top_df, aes(x = EnrichmentScore, y = reorder(Description, EnrichmentScore), fill = Contrast)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.75) +
    theme_minimal(base_size = 12) +
    labs(title = title, x = "-log10(adj p-value)", y = "Term/Pathway", fill = "Contrast")
  
  ggsave(outfile, plot = p, width = 13, height = 10, dpi = 300)
  invisible(p)
}

plot_ora_bar(ora_go,
             "GO BP ORA (top terms per contrast)",
             file.path(out_dir_ora, "ORA_GO_Combined_Barplot.png"),
             top_n = 10)

plot_ora_bar(ora_reactome,
             "Reactome ORA (top pathways per contrast)",
             file.path(out_dir_ora, "ORA_Reactome_Combined_Barplot.png"),
             top_n = 10)

cat("✅ ORA enrichment exported.\n")


############################
# 7) GSEA ENRICHMENT: GO BP + Reactome (Up/Down via NES)
############################
ranked_list <- list()

for (contrast in names(contrast_tables)) {
  logfc_col <- contrast
  padj_col  <- get_padj_col(contrast, colnames(df), offset = 6)
  if (is.na(padj_col)) next
  
  d0 <- df %>%
    transmute(
      Protein = .data[[protein_col]],
      logFC   = suppressWarnings(as.numeric(.data[[logfc_col]])),
      padj    = suppressWarnings(as.numeric(.data[[padj_col]]))
    ) %>%
    filter(!is.na(logFC) & !is.na(padj)) %>%
    # relaxed prefilter to build the ranked list (helps “weak” contrasts)
    filter(abs(logFC) > gsea_prefilter_logfc | padj < gsea_prefilter_padj)
  
  # UniProt -> Entrez
  d0 <- bind_cols(d0, parse_protein_fields(d0$Protein))
  
  entrez <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = d0$UniProtID,
    column = "ENTREZID",
    keytype = "UNIPROT",
    multiVals = "first"
  )
  d0$Entrez <- entrez
  d0 <- d0 %>% filter(!is.na(Entrez))
  
  # remove duplicate Entrez (required: names(stats) must be unique)
  d0 <- d0[!duplicated(d0$Entrez), ]
  
  stats <- d0$logFC
  names(stats) <- d0$Entrez
  stats <- sort(stats, decreasing = TRUE)
  
  ranked_list[[contrast]] <- stats
}

all_gsea_go <- list()
all_gsea_reactome <- list()

for (contrast in names(ranked_list)) {
  cat("Running GSEA for:", contrast, "\n")
  gene_vector <- ranked_list[[contrast]]
  
  if (length(gene_vector) < gsea_min_ranked_genes) {
    message("⚠️ Skipping (too few ranked genes): ", contrast)
    next
  }
  
  gsea_go <- tryCatch({
    clusterProfiler::gseGO(
      geneList     = gene_vector,
      OrgDb        = org.Mm.eg.db,
      keyType      = "ENTREZID",
      ont          = "BP",
      minGSSize    = gsea_min_gs_size,
      maxGSSize    = gsea_max_gs_size,
      pvalueCutoff = gsea_pvalue_cutoff,
      verbose      = FALSE
    )
  }, error = function(e) NULL)
  
  if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
    go_df <- as.data.frame(gsea_go)
    go_df$Contrast <- contrast
    all_gsea_go[[contrast]] <- go_df
  } else {
    message("⚠️ No GO GSEA for: ", contrast)
  }
  
  gsea_re <- tryCatch({
    ReactomePA::gsePathway(
      geneList     = gene_vector,
      organism     = "mouse",
      pvalueCutoff = gsea_pvalue_cutoff,
      minGSSize    = gsea_min_gs_size,
      verbose      = FALSE
    )
  }, error = function(e) NULL)
  
  if (!is.null(gsea_re) && nrow(as.data.frame(gsea_re)) > 0) {
    re_df <- as.data.frame(gsea_re)
    re_df$Contrast <- contrast
    all_gsea_reactome[[contrast]] <- re_df
  } else {
    message("⚠️ No Reactome GSEA for: ", contrast)
  }
}

gsea_go <- bind_rows(all_gsea_go)
gsea_reactome <- bind_rows(all_gsea_reactome)

write.csv(gsea_go,       file.path(out_dir_gsea, "GSEA_GO_All_Contrasts.csv"), row.names = FALSE)
write.csv(gsea_reactome, file.path(out_dir_gsea, "GSEA_Reactome_All_Contrasts.csv"), row.names = FALSE)

# Per-contrast facet plots (top 10 by |NES|)
plot_gsea_facet <- function(df, title, outfile, top_n = 10) {
  if (nrow(df) == 0) return(invisible(NULL))
  
  top_df <- df %>%
    group_by(Contrast) %>%
    slice_max(order_by = abs(NES), n = top_n, with_ties = FALSE) %>%
    ungroup()
  
  p <- ggplot(top_df, aes(x = reorder(Description, NES), y = NES, fill = NES > 0)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~ Contrast, scales = "free_y") +
    coord_flip() +
    scale_fill_manual(values = c("blue", "red")) +
    theme_minimal(base_size = 12) +
    labs(title = title, x = "Gene set", y = "NES (red = up, blue = down)")
  
  ggsave(outfile, plot = p, width = 12, height = 10, dpi = 300)
  invisible(p)
}

plot_gsea_facet(gsea_go,
                "GSEA GO BP (top 10 per contrast)",
                file.path(out_dir_gsea, "GSEA_GO_Facet.png"),
                top_n = 10)

plot_gsea_facet(gsea_reactome,
                "GSEA Reactome (top 10 per contrast)",
                file.path(out_dir_gsea, "GSEA_Reactome_Facet.png"),
                top_n = 10)

# “Single plot” like your reference: pick global top pathways and plot across all contrasts
plot_gsea_single <- function(df, title, outfile, top_global = 40) {
  if (nrow(df) == 0) return(invisible(NULL))
  
  top_terms <- df %>%
    group_by(Description) %>%
    slice_max(order_by = abs(NES), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    slice_max(order_by = abs(NES), n = top_global, with_ties = FALSE) %>%
    pull(Description)
  
  plot_df <- df %>%
    filter(Description %in% top_terms) %>%
    mutate(GeneSet = Description)
  
  p <- ggplot(plot_df, aes(x = NES, y = reorder(GeneSet, NES), fill = Contrast)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.75) +
    theme_minimal(base_size = 13) +
    labs(title = title, x = "NES (positive = up, negative = down)", y = "Gene set", fill = "Contrast") +
    geom_vline(xintercept = 0, linetype = "dashed")
  
  ggsave(outfile, plot = p, width = 16, height = 12, dpi = 300)
  invisible(p)
}

plot_gsea_single(gsea_go,
                 "GSEA GO BP (global top gene sets across contrasts)",
                 file.path(out_dir_gsea, "GSEA_GO_SinglePlot.png"),
                 top_global = 40)

plot_gsea_single(gsea_reactome,
                 "GSEA Reactome (global top pathways across contrasts)",
                 file.path(out_dir_gsea, "GSEA_Reactome_SinglePlot.png"),
                 top_global = 40)

cat("✅ GSEA exported.\n")


############################
# 8) DEG OVERLAP: UpSet + counts
############################
# Use UniProtName lists (or switch to UniProtID)
deg_sets <- lapply(names(contrast_tables), function(contrast) {
  contrast_tables[[contrast]] %>%
    filter(Regulation != "Not significant") %>%
    pull(UniProtName) %>%
    unique()
})
names(deg_sets) <- names(contrast_tables)

if (length(deg_sets) >= 2) {
  pdf(file.path(out_dir_overlap, "DEG_Overlap_UpSet.pdf"), width = 10, height = 6)
  UpSetR::upset(UpSetR::fromList(deg_sets), nsets = length(deg_sets), order.by = "freq")
  dev.off()
  
  deg_summary <- data.frame(
    Contrast = names(deg_sets),
    DEG_Count = sapply(deg_sets, length)
  )
  write.csv(deg_summary, file.path(out_dir_overlap, "DEG_summary_counts.csv"), row.names = FALSE)
  
  cat("✅ UpSet + counts exported.\n")
} else {
  message("⚠️ Not enough DEG sets to draw UpSet.")
}


cat("\n✅ Pipeline finished.\n")
