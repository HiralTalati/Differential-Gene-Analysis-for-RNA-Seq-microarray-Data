#################################################################
# MASTER PIPELINE: OPTIMIZED LONG-FORMAT & MULTI-EXPORT
# Features: 
# - Normalization Audit & QC (Density/UMAP)
# - Global Unbiased DEG (p < 0.05, |logFC| > 0.2)
# - nuID handling (GSE140829) & Max-t Redundancy Collapsing
# - Two Long-Format Master Files (Full vs. Light)
#################################################################

options(download.file.method.GEOquery = "libcurl")
required_packages <- c("GEOquery", "limma", "umap", "ggplot2", "dplyr", 
                       "tidyr", "tibble", "openxlsx", "org.Hs.eg.db", 
                       "illuminaHumanv4.db", "illuminaHumanv3.db", 
                       "lumi", "lumiHumanIDMapping", "AnnotationDbi")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
    library(pkg, character.only = TRUE)
  }
}

run_comprehensive_pipeline <- function(gse_path = NA, dataset_id, group_col) {
  
  # Create directory with dataset name
  out_dir <- paste0(dataset_id, "_Excel_Results")
  dir.create(out_dir, showWarnings = FALSE)
  
  # --- 01: LOAD & EXPLORE ---
  eset_list <- if(is.na(gse_path)) getGEO(dataset_id, destdir = ".", GSEMatrix = TRUE) else getGEO(filename = gse_path, GSEMatrix = TRUE)
  eset <- if(is.list(eset_list)) eset_list[[1]] else eset_list
  expr_mat <- exprs(eset); pheno_df <- pData(eset)
  
  # Clean Groups
  raw_vals <- as.character(pheno_df[[group_col]])
  pheno_df$Group <- case_when(
    raw_vals %in% c("Control", "CTL", "control") ~ "CTL",
    raw_vals %in% c("MCI", "borderline MCI") ~ "MCI",
    raw_vals %in% c("AD", "Alzheimer's Disease") ~ "AD",
    TRUE ~ NA_character_
  )
  pheno_df <- pheno_df[!is.na(pheno_df$Group), ]
  expr_mat <- expr_mat[, rownames(pheno_df)]
  pheno_df$Group <- factor(pheno_df$Group, levels = c("CTL", "MCI", "AD"))
  
  # Print Stats
  cat("\n--- Dataset Exploration:", dataset_id, "---")
  cat("\nTotal Samples:", nrow(pheno_df))
  print(table(pheno_df$Group))
  
  # --- 02: AUDIT & NORMALIZATION ---
  if (max(expr_mat, na.rm=T) > 100) expr_mat <- log2(expr_mat + 1)
  if (max(apply(expr_mat, 2, median)) - min(apply(expr_mat, 2, median)) > 0.5) {
    expr_mat <- limma::normalizeBetweenArrays(expr_mat, method="quantile")
  }
  
  # --- 03: QC PLOTS ---
  # 05: Density
  expr_long_qc <- as.data.frame(expr_mat[1:500,]) %>% rownames_to_column("ID") %>% pivot_longer(-ID)
  ggsave(file.path(out_dir, "05_QC_Density.png"), ggplot(expr_long_qc, aes(x=value, group=name)) + geom_density(alpha=0.1) + theme_minimal())
  
  # 06: UMAP
  set.seed(42); top_p <- order(apply(expr_mat, 1, var), decreasing = T)[1:min(1000, nrow(expr_mat))]
  u_out <- umap::umap(t(expr_mat[top_p, ])); u_df <- data.frame(UMAP1=u_out$layout[,1], UMAP2=u_out$layout[,2], Group=pheno_df$Group)
  ggsave(file.path(out_dir, "06_QC_UMAP.png"), ggplot(u_df, aes(UMAP1, UMAP2, color=Group)) + geom_point(size=2) + theme_minimal())
  
  # --- 04: GLOBAL DEG (SINGLE ACTIVE SHEET LOGIC) ---
  design <- model.matrix(~ 0 + Group, data = pheno_df)
  colnames(design) <- levels(pheno_df$Group)
  fit <- lmFit(expr_mat, design)
  cont_matrix <- makeContrasts(CTL_vs_MCI=MCI-CTL, MCI_vs_AD=AD-MCI, CTL_vs_AD=AD-CTL, levels=design)
  fit2 <- contrasts.fit(fit, cont_matrix)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  
  # --- 05: ANNOTATION & REDUNDANCY (Max-t) ---
  all_probes <- as.character(rownames(expr_mat))
  if (dataset_id == "GSE140829") {
    conn <- lumiHumanIDMapping_dbconn()
    map_tab <- as.data.frame(dplyr::collect(dplyr::tbl(conn, "nuID_MappingInfo"))) %>% 
      dplyr::select(ProbeID = nuID, EntrezID, Symbol) 
  } else {
    lib <- ifelse(dataset_id == "GSE63060", "illuminaHumanv3.db", "illuminaHumanv4.db")
    map_tab <- data.frame(ProbeID = all_probes,
                          EntrezID = as.character(AnnotationDbi::mapIds(get(lib), keys=all_probes, column="ENTREZID", keytype="PROBEID", multiVals="first")),
                          Symbol = as.character(AnnotationDbi::mapIds(get(lib), keys=all_probes, column="SYMBOL", keytype="PROBEID", multiVals="first")))
  }
  map_tab <- map_tab %>% filter(!is.na(EntrezID))
  
  best_probes <- topTable(fit2, coef="CTL_vs_MCI", number=Inf) %>% rownames_to_column("ProbeID") %>%
    inner_join(map_tab, by="ProbeID") %>% group_by(Symbol) %>% slice_max(abs(t), n=1, with_ties=F) %>% pull(ProbeID)
  
  # --- 06: DATA CONSOLIDATION ---
  prec_limit <- 2.22e-16
  master_stats <- data.frame()
  for(cp in colnames(cont_matrix)) {
    tab <- topTable(fit2, coef=cp, number=Inf) %>% rownames_to_column("ProbeID") %>%
      mutate(P.Value = pmax(P.Value, prec_limit), adj.P.Val = pmax(adj.P.Val, prec_limit), Comparison = cp)
    master_stats <- rbind(master_stats, tab)
  }
  
  expr_long <- as.data.frame(expr_mat[best_probes, ]) %>% rownames_to_column("ProbeID") %>%
    pivot_longer(-ProbeID, names_to = "Sample_ID", values_to = "Expression") %>%
    left_join(pheno_df %>% rownames_to_column("Sample_ID") %>% dplyr::select(Sample_ID, Group), by="Sample_ID")
  
  # MASTER FULL TABLE
  final_full <- expr_long %>% inner_join(map_tab, by="ProbeID") %>% left_join(master_stats, by="ProbeID") %>%
    mutate(Significance = case_when(adj.P.Val < 0.05 & logFC >= 0.2 ~ "Upregulated", adj.P.Val < 0.05 & logFC <= -0.2 ~ "Downregulated", TRUE ~ "Not Significant"))
  
  # MASTER LIGHT TABLE (Requested Columns)
  final_light <- final_full %>% dplyr::select(Sample_ID, Group, Symbol, logFC, t, P.Value, adj.P.Val, B)
  
  # --- 07: EXPORTS ---
  write.xlsx(master_stats, file.path(out_dir, "01_Raw_DEG_Stats.xlsx"))
  write.xlsx(final_full, file.path(out_dir, "02_DEG_Annotated_Full_Master.xlsx"))
  write.xlsx(final_light, file.path(out_dir, "03_DEG_Annotated_Light_Master.xlsx"))
  write.xlsx(final_full %>% filter(Significance != "Not Significant"), file.path(out_dir, "04_Significant_DEG.xlsx"))
  
  # Wide Matrix 
  write.xlsx(as.data.frame(expr_mat[best_probes, ]) %>% rownames_to_column("ProbeID") %>% 
               inner_join(map_tab, by="ProbeID") %>% dplyr::select(EntrezID, everything(), -ProbeID, -Symbol), 
             file.path(out_dir, "08_Expression_Matrix_Wide.xlsx"))
  
  # Volcano Plots
  for(cp in unique(master_stats$Comparison)) {
    p_vol <- ggplot(final_full %>% filter(Comparison == cp) %>% distinct(Symbol, .keep_all=T), 
                    aes(x=logFC, y=-log10(adj.P.Val), color=Significance)) + geom_point(alpha=0.4) + 
      scale_color_manual(values=c("Upregulated"="#E41A1C", "Downregulated"="#377EB8", "Not Significant"="#999999")) + 
      theme_minimal() + labs(title=cp)
    ggsave(file.path(out_dir, paste0("09_Volcano_", cp, ".png")), p_vol)
  }
  
  cat("\n[SUCCESS] Completed:", dataset_id)
}

# BATCH RUN
datasets <- list(
  list(id = "GSE63060",  path = "GSE63060_series_matrix.txt.gz", col = "status:ch1"),
  list(id = "GSE63061",  path = "GSE63061_series_matrix.txt.gz", col = "status:ch1"),
  list(id = "GSE140829", path = "GSE140829_series_matrix.txt.gz", col = "diagnosis:ch1") 
)
for (ds in datasets) { tryCatch({ run_comprehensive_pipeline(ds$path, ds$id, ds$col) }, error = function(e) { print(e) }) }
