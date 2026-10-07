
emptyDrops_wrapr <- function(rna_mat = rna_mat, 
                             object = object,
                             lower = 100,
                             niters = 10000) {
  #use a local variable to save the original matrix in here
  original_mat <- rna_mat
  
  initial_cell_cnt <- dim(original_mat)[1]
  initial_gene_cnt <- dim(original_mat)[2]
  
  #start of tryCatch
  tryCatch(
    {
      set.seed(100)
      #lower: Droplets with <100 UMI assumed empty. Depending on the experiment this could change
      #niters: Iterations for specifying the number of iterations to use for the Monte Carlo p-value calculations.
      empty_report <- emptyDrops(
        m = rna_mat,
        lower = 100,
        niters = 10000,
        test.ambient = TRUE
      )
      
      #filtering results
      #Identify cell with FDR < 0.01
      is_cell <- empty_report$FDR < 0.01
      is_cell[is.na(is_cell)] <- FALSE # Treat NAs as empty (very low UMIs per drolet)
      validated_barcodes <- colnames(object)[is_cell] # filtering by validated cells
      
      
      cat("Droplets tested:", ncol(object), "\n")
      cat("Cells called:", sum(is_cell), 
          "(", round(sum(is_cell)/ncol(object)*100, 1), "%)\n")
      cat("Empty droplets removed:", sum(!is_cell), 
          "(", round(sum(!is_cell)/ncol(object)*100, 1), "%)\n")
      
      
      #
      save_raw_counts <- original_mat
      
      #filtering the Seurat Object by the validated BardCodes
      new_object <- subset(object, cells = validated_barcodes)
      
      cat("After Empty Drops:", ncol(new_object), "cells retained\n")
      
      
      return(list(report = empty_report, 
                  subset_object = new_object))
    },
    #tryCatch error function will return the original matrix
    #This part is necessary because the author of the
    error = function(e) {
      message("Warning", e)
      message("Matrix Already Filtered by CellRanger")
      return(original_mat)
    }#end of error bracket
    
  )# end of tryCatch
  
  
}# end of wrapr function


scqc_scatters <- function(tmp_object,
                          dataset_name = "example_dataset",
                          image_title = "pre_threshold"){
  
  qc_images <- "qc_images"
  if(!dir.exists(qc_images)){
    dir.create("qc_images")
  }
  
  QC_testing <- tmp_object
  
  #percentage of mitochondrial RNA counts
  #calculate percentageof mitochondrial RNA per barcode
  QC_testing[["percent.mt"]] <- PercentageFeatureSet(QC_testing, pattern = "^MT-")
  #percentage of ribosomal RNA counts
  QC_testing[["percent.ribo"]] <- PercentageFeatureSet(QC_testing, pattern = "^RP[SL]")
  #What is the complexity of each cell library
  QC_testing$log10GenesPerUMI <- log10(QC_testing$nFeature_RNA) / log10(QC_testing$nCount_RNA)
  
  # Summary statistics to guide threshold setting
  cat("\nQC Metric Distributions:\n")
  cat("nCount_RNA (UMI):\n")
  cat("  Median:", median(QC_testing$nCount_RNA), 
      "| Q1-Q3:", quantile(QC_testing$nCount_RNA, 0.25), "-", 
      quantile(QC_testing$nCount_RNA, 0.75), "\n")
  
  cat("nFeature_RNA (genes):\n")
  cat("  Median:", median(QC_testing$nFeature_RNA),
      "| Q1-Q3:", quantile(QC_testing$nFeature_RNA, 0.25), "-",
      quantile(QC_testing$nFeature_RNA, 0.75), "\n")
  
  cat("percent.mt:\n")
  cat("  Median:", round(median(QC_testing$percent.mt), 2), "%",
      "| 95th percentile:", round(quantile(QC_testing$percent.mt, 0.95), 2), "%\n")
  
  cat("percent.ribo:\n")
  cat("  Median:", round(median(QC_testing$percent.ribo), 2), "%\n")
  
  #setting Idents to the original identity of the object
  Idents(QC_testing) <- "orig.ident"
  
  p4 <- VlnPlot(
    QC_testing,
    features =c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
    ncol = 4,
    pt.size = 0.1
  ) + theme(plot.title = element_text(face = "bold"))
  
  
  #file number 1 for violin plots
  file_1 <-  glue::glue("./qc_images/{dataset_name}_UMI_Feature_MT_scatters_{image_title}.png")
  gg_patchwork(p4, file_1, width = 10, height = 6)
  
  
  #FeatureScatter, UMI vs Gene Detected
  #
  p5 <- FeatureScatter(QC_testing, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    labs(title = "UMI vs Genes Detected")
  
  #FeatureScatter, UMI vs Mitochondrial %
  p6 <- FeatureScatter(QC_testing, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    labs(title = "UMI vs Mitochondrial %")
  
  # Mitochondrial % vs Ribosomal %
  #How many dying cells vs active cells
  p7 <- FeatureScatter(QC_testing, feature1 = "percent.mt", feature2 = "percent.ribo") +
    labs(title = "Mitochondrial % vs Ribosomal %")
  
  p_scatter <- p5 + p6 + p7
  
  file_2 <- glue::glue("./qc_images/{dataset_name}_relation_scatters_{image_title}.png")
  gg_patchwork(p_scatter, file_2, width = 12, height = 6)
  
  message("Visual Credits go Dr. Li Guo of www.ngs101.com")
  
  return(list(vlnplots = p4, 
              rltn_scatters = p_scatter))
}#end of function


scqc_thresholds <- function(tmp_object,
                            dataset_name = "example",
                            nfeature_min = 200,
                            nfeature_max = 7500,
                            ncount_min = 500,
                            mt_thresh = 10){
  
  qc_df <- data.frame(
    nCount_RNA = QC_testing$nCount_RNA,
    nFeature_RNA = QC_testing$nFeature_RNA,
    percent.mt = QC_testing$percent.mt
  )#setting dataframe for qc filters
  
  qc_df$pass_qc <- (
    qc_df$nCount_RNA >= ncount_min &
      qc_df$nCount_RNA <= ncount_max &
      qc_df$nFeature_RNA >= nfeature_min &
      qc_df$nFeature_RNA <= nfeature_max &
      qc_df$percent.mt < mt_thresh
  )# setting boundary filters for qc data
  
  #Threshol
  p8 <- ggplot(qc_df, aes(x = log10(nCount_RNA + 1), y = log10(nFeature_RNA + 1), 
                          color = pass_qc)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_vline(xintercept = log10(c(ncount_min, ncount_max)), 
               linetype = "dashed", color = "red") +
    geom_hline(yintercept = log10(c(nfeature_min, nfeature_max)), 
               linetype = "dashed", color = "red") +
    scale_color_manual(values = c("TRUE" = "#06D6A0", "FALSE" = "#EF476F")) +
    
    labs(
      title = "Cell Filtering Thresholds",
      subtitle = paste0(
        sum(qc_df$pass_qc), " cells pass QC (",
        round(sum(qc_df$pass_qc) / nrow(qc_df) * 100, 1), "%). ",
        "If thresholds are applied, ",
        nrow(qc_df) - sum(qc_df$pass_qc),
        " cells (",
        round((1 - sum(qc_df$pass_qc) / nrow(qc_df)) * 100, 1),
        "%) will be erased."
      ),
      
      x = "log10(UMI + 1)",
      y = "log10(Genes + 1)",
      color = "Pass QC") +  theme_classic()
  
  cat("Visual Credits go Dr. Li Guo of www.ngs101.com")
  
  threshold_image_path <- glue::glue("qc_images/{dataset_name}_threshold_checking.png")
  
  gg_patchwork(p8, threshold_image_path)
  
  return(list(image = p8,
         qc_df = qc_df,
         qc_df_pass = qc_df$pass_qc))
}#end of function





gg_patchwork <- function(plot, filename, width = 8, height = 6, dpi = 300, ...) {
  if (!grepl("\\.png$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".png")
  }
  # Open a device and print the plot explicitly to bypass RStudio window constraints
  grDevices::png(filename, width = width, height = height, units = "in", res = dpi)
  print(plot)  # works for ggplot OR patchwork
  dev.off()
  message("Saved: ", normalizePath(filename))
}