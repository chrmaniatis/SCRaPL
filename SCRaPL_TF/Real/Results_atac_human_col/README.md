Folder used to store inference results on PBMC data and perform integration . 

To replicate integration results use the following steps:

1)Run Human_estimate_cor to produce the correlation matrix between all features.

2)Run Pick_features.ipynd to pick 60k features with the highest in magnitude Pearson correlation.

3)Run Human_make_data to load pre-process and save PBMC data.

4)Run SCRaPL_atac_human(in parent direcotry) to sample parameters of interest.

5)Run Col_prts to aggregate posterior samples form different chunks.

6)Sample latent space from the model using Gen_dat.

7)Run Seurat_my_data_int to perform integration on SCRaPL denoised data.

8)Run Seurat_hum to perform integration on raw data.
