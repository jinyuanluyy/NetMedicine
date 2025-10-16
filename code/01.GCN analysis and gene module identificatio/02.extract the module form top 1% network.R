#co-expression network----------random walk
##extract the module form top 1% network, neworkson-----sunjie code
setwd("./Codes/")
source('networkson_change.R')

top_quantile=0.99 
num_nodes=20 
transitivity=0.6 
exp_th=2

load("CPTAC_PDAC_tumor_exp.Rdata") 
mean_value <- as.matrix(rowMeans(CPTAC_PDAC_tumor_exp))
index <- which(mean_value > exp_th)
symbol_exp <- CPTAC_PDAC_tumor_exp[index,]
path_out <- "result_CPTAC_PDAC_top1"
setwd(path_out)

corMatrix = makeCorTable(symbol_exp, cutoff = top_quantile, mode = "spearman", self=F, debug =F) 
save(file="PDAC_coexp_network_top1.Rdata",corMatrix)

corNet = makeCorNet(corMatrix) 
moduleList = makeModuleList(corNet, debug = F) 
save(file="module_list.Rdata", moduleList)

cytoscapematerial = annotateModulesByCC(corNet, moduleList, cutCluster = num_nodes, cutCC = transitivity, debug = F) # 提取节点数大于阈值的模块

write.node.cytoscape(cytoscapematerial$nodeTable, 'cytoNodefile.txt') 


all_gene_pairs <- data.frame(Gene1 = character(), Gene2 = character(), stringsAsFactors = FALSE)

for (i in seq_along(corNet)) {

  current_element <- corNet[[i]]


  vertex_name <- names(current_element)[1]


  connected_genes <- names(current_element[[1]])


  if (length(connected_genes) > 0) {

    pairs_df <- data.frame(Gene1 = vertex_name, Gene2 = connected_genes, stringsAsFactors = FALSE)


    all_gene_pairs <- rbind(all_gene_pairs, pairs_df)
  }
}



print(all_gene_pairs)


write.table(all_gene_pairs, file = "CPTAC_all_gene_pairs.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


