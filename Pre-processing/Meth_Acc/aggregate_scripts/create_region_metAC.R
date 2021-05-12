############################################################################################################################
##Code is adapted from Melissa: Bayesian clustering and imputation of single-cell methylomes. Author is Andreas Kapourani.##
############################################################################################################################
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(GenomicRanges))

#
# # Data files
io <- list(anno_name = "prom_50_50")#active_enhancers
io$base_dir   <- "~/SCRaPL/Pre-processing/Meth_acc/" #Directory with scripts 
io$out_dir    <- paste0(io$base_dir, "/acc/processed/")#change met to acc for accessibility data
io$annos_file <- paste0(io$base_dir, "/annotations/", io$anno_name, ".bed")
io$met_dir    <- paste0(io$base_dir, "/acc/binarised/")#change met to acc for accessibility data    
io$met_files  <- list.files(io$met_dir, pattern = "*.gz", full.names = TRUE)

#
# # Parameter options
opts <- list()
opts$upstream   <- -5000   # Upstream of region
opts$downstream <- 5000    # Downstream of region
opts$chrom_size <- NULL    # Chromosome size file
opts$cov        <- 3       # Promoters with at least n CpGs
opts$sd_thresh  <- -1      # Threshold for variance of methylation across region

#
# # Read annotation file
annos <- fread(input = io$annos_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
               showProgress = FALSE) %>%
    setnames(c("chr", "start", "end", "strand", "id", "anno")) %>%
    .[, c("anno") := list(NULL)] %>%
    #.[, c("anno", "chr") := list(NULL, as.factor(paste0("chr", chr)))] %>%
    setkey(chr, start, end) #%>% GRanges()

# Create annotation region
anno_region <- copy(annos)
anno_region <- anno_region %>% .[, centre := floor((start + end)/2) ]

# Only for Nanog regions we create a larger genomic region due to sparse CpG coverage
if (io$anno_name != "active_enhancers") {
    anno_region <- anno_region[, c("start", "end") := list(centre + opts$upstream,
                                                           centre + opts$downstream)]
}
# Create GRanges objects
anno_region <- GRanges(anno_region)
annos <- GRanges(annos)

# # Create methylation regions
met <- list()
for (m_file in io$met_files) {
    # Extract cell ID
    # TODO: Should check for other files, since names might change!!!!!!!!!!!!!!!!!!
    cell <-gsub("/home/christos/Desktop/MTDATA/acc//binarised/","\\1",m_file)
    cell <-gsub(".tsv.gz","\\1",cell)
   # cell <- gsub(".*_([^;]+)\\.tsv.*", "\\1", m_file)
    print(cell)
    # Read scBS seq data
    met_dt <- BPRMeth::read_met(file = m_file, type = "sc_seq") %>% as.data.table %>%
      .[, "met" := ifelse(met > .5, 1, 0)] %>% setkey(seqnames, start, end) %>% GRanges()
    # Create promoter methylation regions
    met[[cell]] <- BPRMeth::create_region_object(met_dt = met_dt, anno_dt = anno_region,
                                        cov = opts$cov, sd_thresh = opts$sd_thresh,
                                        ignore_strand = TRUE, filter_empty_region = FALSE)$met
}
rm(met_dt)

#
# # Store the results
message("Storing results...")
obj <- list(met = met, anno_region = anno_region, annos = annos,
            io = io, opts = opts)
saveRDS(obj, file = paste0(io$out_dir, "unfiltered/", io$anno_name, ".rds"))
  