# plot/table of each cohort + cancer_group or cancer_group
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

tumor_plot <- function(expr_mat_gene, hist_file, 
                       analysis_type = c("cohort_cancer_group_level", "cancer_group_level"), 
                       log = FALSE, 
                       plots_dir = 'www', results_dir = 'www', 
                       plot_width = 10, plot_height = 9, mapping_file = 'www/metadata.tsv'){
  
  # subset histology to minimal columns
  hist_file <- hist_file %>%
    select(Kids_First_Biospecimen_ID, cohort, sample_type, cancer_group)
  
  # expression matrix to long format
  expr_mat_gene <- expr_mat_gene %>%
    tidyr::gather("Kids_First_Biospecimen_ID", "tpm", -c("gene"))
  
  # combine with histology file
  expr_mat_gene <- expr_mat_gene %>%
    inner_join(hist_file, by = "Kids_First_Biospecimen_ID")
  
  # format x-axis labels and filter to n >= 5
  if(analysis_type == "cohort_cancer_group_level"){
    expr_mat_gene <- expr_mat_gene %>%
      group_by(cohort, cancer_group) %>%
      mutate(n_samples = n()) %>%
      filter(n_samples >= 5) %>%
      mutate(x_labels = paste0(cancer_group, ", ", cohort,  " (N = ", n_samples, ")"))
  } else if(analysis_type == "cancer_group_level") {
    expr_mat_gene <- expr_mat_gene %>%
      group_by(cancer_group) %>%
      mutate(n_samples = n()) %>%
      filter(n_samples >= 5) %>%
      mutate(x_labels = paste0(cancer_group, " (N = ", n_samples, ")"))
  }
  
  # reorder by median tpm
  fcts <- sort(unique(expr_mat_gene$x_labels))
  expr_mat_gene$x_labels <- factor(expr_mat_gene$x_labels, levels = fcts)
  
  # create unique title and filenames
  gene_name <- unique(expr_mat_gene$gene)
  tumor_cohort <- paste0(unique(expr_mat_gene$cohort), collapse = ", ")
  tumor_cohort_fname <- paste0(unique(expr_mat_gene$cohort), collapse = "_")
  if(analysis_type == "cohort_cancer_group_level"){
    title <- paste(gene_name, "Gene Expression across cohorts", sep = "\n")
  } else {
    title <- paste(gene_name, "Gene Expression across cancers", sep = "\n")
  }
  fname <- paste(gene_name, tumor_cohort_fname, "pan_cancer", analysis_type, sep = "_")
  plot_fname <- paste0(fname, '.pdf')
  table_fname <- paste0(fname, '.tsv')
  
  # data-frame for mapping output filenames with info
  mapping_df <- data.frame(gene = gene_name, 
                           plot_type = "tumors_only", 
                           cohort = tumor_cohort,
                           cancer_group = NA,
                           analysis_type = analysis_type, 
                           plot_fname = plot_fname,
                           table_fname = table_fname)
  mapping_file <- file.path(results_dir, 'metadata.tsv')
  if(!file.exists(mapping_file)){
    write.table(x = mapping_df, file = mapping_file, sep = "\t", row.names = F, quote = F)
  } else {
    write.table(x = mapping_df, file = mapping_file, sep = "\t", row.names = F, col.names = F, quote = F, append = TRUE)
  }

  # output table of gene, median and sd
  output_table <- expr_mat_gene %>%
    group_by(gene, x_labels) %>%
    summarise(mean = mean(tpm),
              median = median(tpm),
              sd = sqrt(var(tpm))) %>%
    mutate(mean = round(mean, digits = 2),
           median = round(median, digits = 2),
           sd = round(sd, digits = 2))
  write.table(x = output_table, file = file.path(results_dir, table_fname), sep = "\t", row.names = F, quote = F)
  
  # boxplot
  cols <- c("Tumor" = "grey80")
  if(log == TRUE){
    expr_mat_gene <- expr_mat_gene %>%
      mutate(tpm = log2(tpm + 1))
    y_lab <- "log2 (TPM)"
  } else {
    y_lab <- "TPM"
  }
  output_plot <- ggplot(expr_mat_gene, aes(x = x_labels, y = tpm, fill = sample_type)) +
    stat_boxplot(geom ='errorbar', width = 0.2, lwd = 0.3) +
    geom_boxplot(lwd = 0.3, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) +
    ylab(y_lab) + xlab("") +
    theme_Publication2(base_size = 12) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle(title) +
    scale_fill_manual(values = cols) + theme(legend.position='none')
  ggsave(plot = output_plot, filename = file.path(plots_dir, plot_fname), device = "pdf", width = plot_width, height = plot_height)
  
  # now remove error bars because plotly adds it by default
  output_plot$layers[[1]] <- NULL
  
  # add labels to outliers
  expr_mat_gene <- expr_mat_gene %>% 
    group_by(x_labels) %>% 
    mutate(OutlierFlag = ifelse((tpm < quantile(tpm, 1/3, na.rm = T) - 1.5*IQR(tpm, na.rm = T)) | 
                                  (tpm > quantile(tpm, 2/3, na.rm = T) + 1.5*IQR(tpm, na.rm = T)), 'Outlier', 'NotOutlier'))%>%
    group_by()
  
  # add geom point layer
  output_plot <- output_plot + 
    geom_point(data = expr_mat_gene %>% 
                 filter(OutlierFlag == "Outlier"), aes(group = x_labels, label = Kids_First_Biospecimen_ID), stroke = 0.3)
  
  # convert to plotly
  output_plot <- ggplotly(output_plot, tooltip = c("label"))
  for (i in 1:length(output_plot$x$data)){
    if (output_plot$x$data[[i]]$type=="box"){
      output_plot$x$data[[i]]$marker$opacity = 0  
      output_plot$x$data[[i]]$hoverinfo = "none"
    }
  }
  
  return(list(output_plot = output_plot,
              output_plot_fname = plot_fname,
              output_table = output_table))
}
