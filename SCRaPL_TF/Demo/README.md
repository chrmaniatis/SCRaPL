File containing all demo scripts.

Demo_10X: 10X synthetic data generation and inference.
Demo_NMT: NMT synthetic data generation and inference.

Demo_10X_large_data: Similar to Demo_10X but intended for large sets of features. Uses independence across features to reduce the dimensionality of the problem.  

Demo_NMT_large_data: Similar to Demo_NMT but intended for large sets of features. Uses independence across features to reduce the dimensionality of the problem.  

Col_prts: Script used to aggregate posterior samples from different chunks into single files. This should be used after Demo_NMT_large_data or Demo_10X_large_data, before proceeding with further analysis.

Constr_neg_exp: Script used to construct negative control data. For a detailed description, please have a look in the Supplementary Materials section S2. 

Gam_thr_est: Script used to estimate gamma threshold from negative control data. Later gamma is used in equation (9) in the manuscript.

Ft_detection: Script used to detect features with strong correlation with SCRaPL, Pearson and Spearman. 

To get features with strong correlation, file should run in the following order:

1)Demo_NMT_large_data,Demo_NMT,Demo_10X_large_data,Demo_10X --> Estimate and save posterior parameters. If SCRaPL parameters exceed 80k use large to avoid memory or convergence issues.
2)If inferred parameters >=80k use Col_prts to aggregate samples, the create a negative ontrol datset using Constr_neg_exp. Otherwise, proceed with Constr_neg_exp.
3)Demo_NMT_large_data,Demo_NMT,Demo_10X_large_data,Demo_10X --> Estimate and save posterior parameters for negative control data
4)Use Gam_thr_est to estimate gamma, which is used for feature detection.
5)Detect features with strong correlation using Ft_detection. 
