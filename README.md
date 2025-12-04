# Impact of a secondary bile acid and Vitamin D on detoxification-related gene expression within enteric cells in colon cancer ![bsc](https://img.shields.io/badge/Dissertation-Completed-brightgreen) ![status](https://img.shields.io/badge/Code-Stable-blue)
## R scripts supporting a completed investigative dissertation project on bile acid– and vitamin D–regulated detoxification pathways in enteric cells.

This repository contains the full R code used for data preprocessing, statistical analysis, and visualisation in my dissertation research at Ulster University. The study examines how lithocholic acid (LCA), 3-keto-LCA, and vitamin D modulate the expression of detoxification-related genes in intestinal epithelial cell lines relevant to colorectal cancer biology. Specifically, it investigates the activity of the vitamin D receptor (VDR) and pregnane X receptor (PXR)—two nuclear receptors central to bile-acid-dependent detoxification—and their downstream effects on metabolic, transporter, and xenobiotic-processing genes.

Circulating bile acids and their metabolic derivatives act as signalling molecules capable of initiating inflammatory, metabolic, and carcinogenic pathways.
Secondary bile acids such as LCA are potent ligands of nuclear receptors (VDR, PXR, FXR), and dysregulation of these pathways is implicated in colon cancer progression. Vitamin D is hypothesised to counteract these effects by activating VDR-dependent transcriptional networks involved in detoxification, metabolism, and cellular protection.

### Objectives
The study aimed to: quantify VDR- and PXR-mediated transcriptional responses after ligand exposure.
Assess induction of phase I, II, and III detoxification genes.
Compare basal vs ligand-induced gene expression in LS180 and Caco-2 cell lines.
Evaluate temporal dynamics using a single-exposure (24 h) experiment and a 12 h / 24 h / 48 h time-course experiment.
Experimental design
Two RT-qPCR-based experiments were performed:
#### Experiment 1 – Single 24 h exposure
Cells were treated with:
Vitamin D (1α,25-dihydroxyvitamin D₃)
Lithocholic acid (LCA)
3-keto-LCA
Vehicle controls
#### Experiment 2 – Time-course (12, 24, 48 h)
Ligand exposures were repeated across three time points to capture regulatory kinetics.

## Data analysis pipeline
All analyses were performed in R using the tidyverse ecosystem.

Importing and cleaning raw Ct data
Calculating ΔCt and ΔΔCt
Computing fold changes relative to vehicle controls
Statistical testing using paired Student’s t-tests
Visualisation of gene expression changes with significance annotation (* p < 0.05, ** p < 0.01, *** p < 0.001).

### Citation
If you use this code or adapt the workflow, please cite:
Juhász, Á. J. (2022). Impact of a secondary bile acid and vitamin D on detoxification-related gene expression within enteric cells in colon cancer. Ulster University. GitHub repository:
https://github.com/agnjuh/rtpcr_detox_genes_vitD_LCA
