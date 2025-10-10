CPTAC_moduel <- read.table("cytoNodefile.txt", header = TRUE)
CPTAC_moduel <- CPTAC_moduel[CPTAC_moduel$cc > 0.6, ]
load("module_list.Rdata")
protective signature <- read.table("common_protective signature.txt")
pro-oncogenic_signatures <- read.table("common_pro-oncogenic_signatures.txt")


hypergeometric_test <- function(selected_genes, all_genes, module_genes) {
  overlap_genes <- intersect(module_genes, selected_genes)  
  overlap <- length(overlap_genes)
  module_size <- length(module_genes)
  total_genes <- length(all_genes)
  selected_size <- length(selected_genes)

  p_value <- phyper(overlap - 1, module_size, total_genes - module_size, selected_size, lower.tail = FALSE)
  return(list(overlap = overlap, overlap_genes = overlap_genes, p_value = p_value))
}


all_genes <- unique(unlist(moduleList[as.character(CPTAC_moduel$node)]))


results_list <- list()


for (i in seq_along(moduleList[as.character(CPTAC_moduel$node)])) {
  module_name <- names(moduleList[as.character(CPTAC_moduel$node)])[i]
  module_genes <- moduleList[[module_name]]
  if (length(module_genes) > 20) { 
    test_result <- hypergeometric_test(pro-oncogenic_signatures$V1, all_genes, module_genes)
    module_total_genes <- length(module_genes)
    module_id <- paste0("M", length(results_list) + 1) 
    results_list[[module_id]] <- data.frame(
      ID = module_id,
      module = module_name,
      total_genes = module_total_genes,
      overlap = test_result$overlap,
      module_genes = paste(module_genes, collapse = ","),  
      overlap_genes = paste(test_result$overlap_genes, collapse = ","),  
      p_value = test_result$p_value
    )
  }
}

CPTAC_pro-oncogenic_module_p_values <- do.call(rbind, results_list)

results_list <- list()

for (i in seq_along(moduleList[as.character(CPTAC_moduel$node)])) {
  module_name <- names(moduleList[as.character(CPTAC_moduel$node)])[i]
  module_genes <- moduleList[[module_name]]
  if (length(module_genes) > 20) { 
    test_result <- hypergeometric_test(protective signature$V1, all_genes, module_genes)
    module_total_genes <- length(module_genes)
    module_id <- paste0("M", length(results_list) + 1) 
    results_list[[module_id]] <- data.frame(
      ID = module_id,
      module = module_name,
      total_genes = module_total_genes,
      overlap = test_result$overlap,
      module_genes = paste(module_genes, collapse = ","), 
      overlap_genes = paste(test_result$overlap_genes, collapse = ","),  
      p_value = test_result$p_value
    )
  }
}


CPTAC_protective_module_p_values <- do.call(rbind, results_list)


filtered_CPTAC_pro-oncogenic_module_p_values <- CPTAC_pro-oncogenic_module_p_values %>%
  filter(p_value < 0.05)


filtered_CPTAC_pro-oncogenic_module_p_values$Module_Features <- "pro-oncogenic_module"


filtered_CPTAC_protective_module_p_values <- CPTAC_protective_module_p_values %>%
  filter(p_value < 0.05)


filtered_CPTAC_protective_module_p_values$Module_Features <- "protective_module"


filtered_CPTAC_pro-oncogenic_module_p_values_insignificant <- CPTAC_pro-oncogenic_module_p_values %>%
  filter(!ID %in% filtered_CPTAC_pro-oncogenic_module_p_values$ID) %>%
  filter(!ID %in% filtered_CPTAC_protective_module_p_values$ID)


filtered_CPTAC_pro-oncogenic_module_p_values_insignificant$Module_Features <- "insignificant"


CPTAC_module<- bind_rows(filtered_CPTAC_pro-oncogenic_module_p_values, filtered_CPTAC_protective_module_p_values, filtered_CPTAC_pro-oncogenic_module_p_values_insignificant)

CPTAC_module <- CPTAC_module %>%
  mutate(
    overlap = if_else(Module_Features == "insignificant", NA_real_, overlap),
    overlap_genes = if_else(Module_Features == "insignificant", NA_character_, overlap_genes),
    p_value = if_else(Module_Features == "insignificant", NA_real_, p_value)
  )



save(CPTAC_module, file = "CPTAC_module.Rdata")







