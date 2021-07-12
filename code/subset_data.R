# script to reduce file sizes by limiting to current cohorts
# read input data
expr_mat <- readRDS('data/gene-expression-rsem-tpm-collapsed.rds')
hist_file <- data.table::fread('data/histologies.tsv')

# filter histologies 
hist_file <- hist_file %>%
  filter(experimental_strategy == "RNA-Seq",
         sample_type == "Normal" | !is.na(cancer_group),
         Kids_First_Biospecimen_ID %in% colnames(expr_mat),
         cohort %in% c("GMKF", "GTEx", "PBTA"))
write.table(hist_file, file = 'data/histologies.tsv', quote = F, sep = "\t", row.names = F)

# match expression matrix to histologies file
expr_mat <- expr_mat %>%
  select(hist_file$Kids_First_Biospecimen_ID)

# subset to 10 NBL genes
expr_mat <- expr_mat[grep("^ALK$|^MYCN$|^GPC2$|^CAMKV$|^L1CAM$|^NCAM1$|^B4GALNT1$|^DLL3$|^PHOX2B$|^DLK1$", rownames(expr_mat)),]
saveRDS(expr_mat, file = 'data/gene-expression-rsem-tpm-collapsed-subset.rds')
