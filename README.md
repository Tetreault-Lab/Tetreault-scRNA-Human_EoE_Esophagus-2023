# Tetreault-scRNA-Human_EoE_Esophagus-2023

This repository contains the R code for processing the raw data and for the analyses that were generated for the Tetreault Lab paper Suprabasal cells retain progenitor cell identity programs in eosinophilic esophagitis-driven basal cell hyperplasia. JCI Insight. 2023. https://doi.org/10.1172/jci.insight.171765. Raw data and matrices are available in GEO under accession code GSE218607. Analyses were done using Seurat (v4.2.0) and Monocle3 (v1). Additional libraries are specified throughout.


**There are 6 walkthroughs:**

    1. PreProcess_and_Cluster (Initial pre-processing and scRNA-Seq QC using Seurat, as well as integration and sub-clustering of esophageal epithelial cells)

    2. Downstream_Seurat_Analyses (Differential expression analyses and other downstream Seurat analyses)

    3. Monocle_Analyses (Downstream monocle3 pseudotime analysis for merged object, learned graph, and pseudotime-dependent gene module identification)

    4. Plots_for_Figures (generation of plots used in the paper)

    5. Plot_Functions (functions written for data visualization, these functions are used in other walkthroughs)

    6. GERD_analyses (Projection of GERD esophageal epithelial cells)

**Major Requirements:**

Seurat (4.2.0): https://satijalab.org/seurat/

Monocle (3.0.0): https://cole-trapnell-lab.github.io/monocle3/
