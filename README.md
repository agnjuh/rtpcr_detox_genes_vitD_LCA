This repository contains the R code developed as part of my dissertation project at Ulster University (2022). The project investigates the role of vitamin D receptor (VDR) and pregnane X receptor (PXR) in the regulation of detoxification genes within intestinal cell lines in response to endocrine and bile acid-based ligand activation.

Project overview: circulating bile acids (BA) play a key role in initiating inflammatory processes and tumorigenesis. Vitamin D is hypothesised to inhibit these pathways through activation of VDR, which, along with PXR, acts as a sensor for endogenous toxic metabolites. This study aims to:

assess VDR and PXR activity through gene expression profiling.
evaluate the transactivation of metabolic and transporter genes via ligand exposure.
compare basal and induced gene expression in LS180 and Caco-2 intestinal cell lines.
Experimental setup Two independent experiments were conducted using RT-qPCR:

experiment 1: single 24-hour ligand exposure.
experiment 2: time-course analysis at 12, 24, and 48 hours. Target genes included phase I, II, and III metabolic genes, measured across treated and vehicle control conditions.
Data analysis pipeline The analysis was performed in R using Tidyverse packages. Main steps include:

processing raw CT values.
calculating ΔCT and ΔΔCT using the Livak & Schmittgen (2001) method.
visualising fold change in expression relative to vehicle controls.
performing paired Student's t-tests for statistical significance.
annotating significance levels (*p < 0.05, **p < 0.01, ***p < 0.001) on plots.
Citation:

This work was conducted by Ágnes Judit Juhász as part of the BSc Hons dissertation in Biomedical Science at Ulster University. If you use this code or data in your work, please cite:

Juhász, Á. J. (2022). Impact of a secondary bile acid and vitamin D on detoxification related gene expression within enteric cells in colon Cancer [GitHub repository]. Retrieved from https://github.com/agnjuh/rtpcr_detox_genes_vitD_LCA

@misc{juhasz2022detox,
  author       = {Ágnes Judit Juhász},
  title        = {Impact of a secondary bile acid and vitamin D on detoxification related gene expression within enteric cells in colon cancer},
  year         = {2022},
  url          = {https://github.com/agnjuh/rtpcr_detox_genes_vitD_LCA},
  note         = {GitHub repository}
}
