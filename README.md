Single-cell multi-omics assays offer unprecedented opportunities to explore epigenetic regulation at cellular level. However, high levels of technical noise and data sparsity frequently lead to a lack of statistical power in correlative analyses, identifying very few, if any, significant associations between different molecular layers. Here we propose SCRaPL, a novel computational tool that increases power by carefully modelling noise in the experimental systems. We show on real and simulated multi-omics single-cell data sets that SCRaPL achieves higher sensitivity and better robustness in identifying correlations, while maintaining a similar level of false positives as standard analyses based on Pearson and Spearman correlation.

The repository contains a Matlab and a Tensorflow implementation of SCRaPL. Users are encouraged to choose the Tensorflow implementation which is more scalable. 

Pre-processing: Preprocessing scripts. This folder contains R scripts used to aggregate data from cells for genomic regions of interest,to join epigenomic with transcriptomic layers, to perform scRNA-seq normalization and QC.

SCRaPL_TF: Tensorflow implementration of SCRaPL. This folder contains three subfolders. One for synthetic data analysis, one for analysis on real data and a file with demos.

Paper/Sketches: Backup of inkscape files.

Analysis: Matlab scripts used to produce analysis in the results section.

Bayesian_Inference: Matlab scripts used for inference.

Negative_Control: Matlab scripts used to construct negative control datasets.

Synth_analysis: Matlab scripts used to perform synthetic data analysis.
