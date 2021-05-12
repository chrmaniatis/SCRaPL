############################################################################################################################
##Code is adapted from Melissa: Bayesian clustering and imputation of single-cell methylomes. Author is Andreas Kapourani.##
############################################################################################################################
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(stringr))

#data_file  <- "prom_50_50"
#name_file <- "prom_5000_5000"
data_file <- "active_enhancers"
name_file <- "active_enanhcers_1"
data_dir   <- "~/SCRaPL/Pre-processing/Meth_Acc/met/processed/unfiltered/"
out_dir <- "out"

# Load data
obj <- readRDS(paste0(data_dir, data_file, ".rds"))
met <- obj$met
anno_region <- obj$anno_region
annos <- obj$annos

# Update options
opts <- obj$opts
opts$cov <- 3       # CpG density at each source
opts$is_binomial <- TRUE # If we need to return a big vector of all zeros and 1s

# Consider only regions with enough CpG coverage
met <- lapply(met, function(x) lapply(x, function(y){
    if (NROW(y) < opts$cov) return(NA) else return(y) }))

# If not binomial
if (!opts$is_binomial) {
    # dt <- data.table(met_level = numeric(), expr_level = numeric(), cell = numeric())
    dt <- matrix(0, ncol = 3, nrow = NROW(anno_region) * length(met))
    iter <- 1
    # Iterate over each gene
    for (i in 1:NROW(anno_region)) {
        if (i %% 100 == 0) { print(i) }
        # Iterate over each cell
        for (c in 1:length(met)) {
            # Compute methylation level, -1 for missing data
            met_dt <- met[[c]][[i]]
            if (is.na(met_dt)) {
                met_level <- -1
            }else{
                met_level <- mean(met_dt[,2])
            }
            # Extract gene expression level
            #expr_level <- rna[i, c + 1, with = FALSE]
            dt[iter, ] <- c(met_level, anno_region$id[i], names(met)[c])
            # rbind(dt, data.table(met_level = met_level, expr_level = as.numeric(expr_level), cell = c))
            iter <- iter + 1
        }
    }
    dt2 <- data.table(met = dt[,1], feature = dt[,2], cell = dt[,3])

    fwrite(dt2, file = paste0(data_dir, "out/", name_file, ".csv"), sep = "\t")
} else{
    # dt <- data.table(met_level = numeric(), expr_level = numeric(), cell = numeric())
    dt <- matrix(0, ncol = 4, nrow = NROW(anno_region) * length(met))
    iter <- 1
    # Iterate over each gene
    for (i in 1:NROW(anno_region)) {
        if (i %% 100 == 0) { print(i) }
        # Iterate over each cell
        for (c in 1:length(met)) {
            # Compute methylation level, -1 for missing data
            met_dt <- met[[c]][[i]]
            if (is.na(met_dt)) {
                met_cpgs <- -1
                total_cpgs <- -1
            }else{
                met_cpgs <- length(which(met_dt[,2] == 1))
                total_cpgs <- length(met_dt[,2])
            }
            # Extract gene expression level
            #expr_level <- rna[i, c + 1, with = FALSE]
            dt[iter, ] <- c(met_cpgs, total_cpgs, anno_region$id[i], names(met)[c])
            # rbind(dt, data.table(met_level = met_level, expr_level = as.numeric(expr_level), cell = c))
            iter <- iter + 1
        }
    }
    dt2 <- data.table(met_cpgs = dt[,1], total_cpgs = dt[,2], feature = dt[,3], cell = dt[,4])

    fwrite(dt2, file = paste0(data_dir, "out/", name_file, ".csv"), sep = "\t")
}
