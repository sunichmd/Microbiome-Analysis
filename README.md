## Microbiome Analysis Project (English Version)

This project consolidates and optimizes multiple workflows for microbiome data analysis, mainly covering:



### EasyAmplicon

This workflow is modified from the [YongxinLiu/EasyAmplicon](https://github.com/YongxinLiu/EasyAmplicon) repository (version: April 2025).
The original project was released under the GPL v3 license, which is also adopted here. On top of that, this version introduces enhanced visualization and output modules.

#### Updates to Qiime2 Pipeline

* **Database Updates**

  * Added `silva` database (Qiime2 v2024.10 compatible): [Download link](https://figshare.com/ndownloader/files/56053826)
  * Added `greengene2_24_09` database: [Download link](https://figshare.com/ndownloader/files/56053796)
  * Added `HOMD v15.23` database for oral microbiome

* **Pipeline Optimization (pipeline_qiime2_4ana.sh)**

  * Added: `export UNIFRAC_MAX_CPU=basic` to handle AVX2-incompatible CPUs
  * Updated: taxonomy classifier training script
  * Added: conversion of Qiime2 output into human-readable tables
  * Added: R-based visualization scripts for Qiime2 results
  * Added: LEfSe-compatible format conversion + LEfSe analysis
  * Added: longitudinal analysis (alpha/beta diversity, taxa abundance)
  * Added: alpha.txt and beta.txt conversion into `.qza` format

#### Usearch + Vsearch Pipeline (Windows & Linux)

* Refined pipeline into `pipeline_4ana.sh`

* Optional reference DBs:

  * `silva-138-99-seqs.usearch.fa`: [Download](https://figshare.com/ndownloader/files/56056067)
  * `greengene2_2022.10.seqs.usearch.fna`: [Download](https://figshare.com/ndownloader/files/56056076)

* Visualization Scripts:

  * `alpha_boxplot1.R`: customizable alpha diversity plots with grouping
  * `beta_pcoa2.R`: PCoA + boxplot with group-wise separation tests
  * `tax_stackplot1.R`: customizable stacked plots with faceting
  * `tax_circlize_xlim90.R`: circular taxonomic plots with rotated labels
  * `compare_volcano1.R`: improved volcano plots with correct coloring
  * `compare_heatmap1.sh`: heatmap with optional row annotation
  * **Appendix 2.1**: support for `.fasta` input sequence files

#### Vsearch Pipeline (Mac OS)

* Pipeline adjustments:

  * Replace `cat -A` with `cat -vet` for macOS compatibility
  * Added: R-based rarefaction curve (alpha diversity)
  * Added: R-based beta diversity analysis


### Advanced Analysis Modules

* `BatchCorrect`: Evaluate and correct batch effects using MMUPHin, ComBat, removeBatchEffect, Percentile Normalization
* `SourceTrackerFeast`: Bug-fixed SourceTracker1.r, added FEAST1.Rmd for individual-level source tracking over time
* `Multi-factorAnalysis`: Alpha/beta/taxa comparisons under multifactor designs, supports Qiime2 longitudinal mode
* `RandomForestClassification`: RF classification with group-wise sampling, AUC curve, full-sample heatmap
* `Net_SingleMulti_Vis`: One-command correlation network analysis and visualization with qgraph/Cytoscape export

  * Single-file or two-file correlation network construction
  * Group-wise correlation network calculation
  * Network property metrics (e.g., edge number, density)
  * ZiPi index calculation
  * Module (community) detection
  * Network stability evaluation
  * Layouts using qgraph-based pseudo-graph structure
  * Export-ready results for Graphi visualization
  * Export-ready results for Cytoscape visualization
  * R-based automated Cytoscape plotting
  * Core taxa identification and statistics
* `sparcc_vis_as_cellPaper.Rmd`: Implements SparCC-based correlation network construction and visualization in R, following the layout and styling used in Cell publications.
* `cor_2file.Rmd`: Computes basic correlation networks between features across two datasets, and generates interactive visual outputs.
* `Phytopathogen_Bac_predict`: Pathogenic bacteria prediction + differential analysis
* `calculate_core_otu_by_usearch_treeplot.Rmd`: Core OTU detection + phylogenetic tree
* `adonis_lost_variable_solution.Rmd`: Fix variable-dropping issues in Adonis models
* `calculate_each_tax_raw_count.Rmd`: Compute alpha/beta diversity at taxonomic levels
* `PCoA_density_compare_multi-factor.Rmd`: PCoA density + statistical comparisons
* `Fungi_function_predict.Rmd`: Fungal functional annotation
* `calculate_Bray-Curtis_dissimilarities.Rmd`: Bray-Curtis computation
* `Other_ana.Rmd`: Miscellaneous tools

  * Adonis models
  * LEfSe in R
  * Ternary plots (for 3-group data)
  * NMDS ordination
  * Format greengene2 for usearch
  * Relative abundance-based diversity
  * Taxa accumulation curve
  * Relative abundance impact plots
  * SVD analysis

---

### EasyMetagenome

Modified from [YongxinLiu/EasyMetagenome](https://github.com/YongxinLiu/EasyMetagenome) (version: Dec 2024), this version improves file structure and compatibility.


### EasyMicrobiome

Based on [YongxinLiu/EasyMicrobiome](https://github.com/YongxinLiu/EasyMicrobiome) (version: Dec 2024), this version has been partially refactored to support batch processing.

#### Key Updates:

* Added greengene2 (2022.10) and silva-138 databases into `usearch` folder
* Added `alpha_calculate_pd.R` script: computes 7 types of alpha diversity (requires phylogenetic tree for PD index)

---

> All derivative workflows retain original authorship and licenses. Modifications are clearly documented with the goal of improving reusability and modularity in microbiome analysis.

---