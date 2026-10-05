
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
      
      
      return(list(empty_droplets, new_object))
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