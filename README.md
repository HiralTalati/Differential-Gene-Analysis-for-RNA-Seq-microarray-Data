**Transcriptomic Master Pipeline: Optimized Long-Format & Multi-Export**

This repository contains an end-to-end bioinformatics pipeline written in R for downloading, preprocessing, analyzing, and structuring microarray transcriptomic data. The pipeline is optimized to handle Illumina-based clinical cohorts (specifically targeting Alzheimer's Disease, Mild Cognitive Impairment, and Control setups) using multi-tier normalization audits, probe redundancy collapse, and advanced long-format data reshaping.

📌 Features

• Multi-Platform Dependency Automation: Dynamically audits, downloads, and mounts required CRAN and Bioconductor packages (limma, GEOquery, lumi, database mapping wrappers).

• Normalization Audit & QC: Detects cross-sample scale variances, automatically applies log₂ scaling and quantile normalization between arrays when necessary, and generates UMAP and sample distribution density tracking plots.

• Illumina nuID & Probe Mapping: Resolves complex Illumina ID structures, supporting standard probe keys (v3/v4 architectures) as well as advanced nuID strings using SQLite databases.

• Max-t Redundancy Collapsing: Discards multi-mapping array artifact probes by selecting the high-confidence probe showing the maximum absolute t-statistic per unique Gene Symbol.

• Dual Long-Format Generation: Reshapes millions of wide matrix expression data metrics into highly querying-friendly flat layouts (Full Data Model vs. Light Downstream Model).

🛠️ Tech Stack & Dependencies

The computational workflow runs natively inside an R environment (v4.0+).
Bioconductor & CRAN Packages:
R

# Managed automatically inside the script execution cycle
GEOquery, limma, umap, ggplot2, dplyr, tidyr, tibble, openxlsx, 
org.Hs.eg.db, illuminaHumanv4.db, illuminaHumanv3.db, lumi, 
lumiHumanIDMapping, AnnotationDbi

📂 Input Configurations & Batch Layout
The pipeline executes sequentially over specified targeted NCBI Gene Expression Omnibus (GEO) matrix archives. Ensure the local target .txt.gz family matrices are placed within your root workspace path:

🚀 Execution Workflow

1. Group Filtering & Harmonization
Raw categorical patient clinical arrays are cleaned and standardized into three hard factors: CTL, MCI, and AD. Extraneous labels or missing entries are safely filtered away before computing contrast configurations.

2. Differential Expression Analysis (DEG)
Linear modeling is executed via limma::lmFit across structured matrix contrasts:
• CTL vs MCI (Mild Cognitive Impairment Progression)
• MCI vs AD (Late-stage Transition)
• CTL vs AD (Global Pathological Alterations)
Empirical Bayes (eBayes) statistics are run with structural trends and robust standard error modifications enabled.

3. Pipeline Run Command
Open an R session or terminal interface and execute the batch loop:
bash
Rscript master_pipeline_v21.R

📊 Exported Workspace Assets
For every dataset processed, a separate results folder named {Dataset_ID}_Excel_Results/ is generated, containing the following validation plots and spreadsheets:
text
├── GSE63060_Excel_Results/
│   ├── 05_QC_Density.png               # Diagnostic sample curve distributions
│   ├── 06_QC_UMAP.png                  # Unbiased spatial group clustering 
│   ├── 01_Raw_DEG_Stats.xlsx           # Complete uncollapsed contrast matrices
│   ├── 02_DEG_Annotated_Full_Master.xlsx# Long-format combined metadata (Expression + Stats)
│   ├── 03_DEG_Annotated_Light_Master.xlsx# Condensed long-format layout for clean modeling
│   ├── 04_Significant_DEG.xlsx         # Extracted boundaries (adj.p < 0.05 & |logFC| > 0.2)
│   ├── 08_Expression_Matrix_Wide.xlsx  # Wide numerical frame tailored for ML / XGBoost input
│   ├── 09_Volcano_CTL_vs_MCI.png       # Standard logFC vs Log10 p-value distributions
│   ├── 09_Volcano_MCI_vs_AD.png
│   └── 09_Volcano_CTL_vs_AD.png

📄 License
This analysis framework is open-source and distributed under the MIT License.
