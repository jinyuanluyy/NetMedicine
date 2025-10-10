library(survival)

path_input<-"./Example_data/"
setwd(path_input)

cancerType<-"CPTAC_PDAC"
exp<-as.matrix(read.table(paste0(cancerType,"_trans_exp_TPM_1.txt"),sep = '\t'))

exp<-exp[,-c(1:7)]
gene_list<-as.matrix(colnames(exp))


dataDir=path_input
path_out_pdf<-"output_pdf_TCGA"
path_out_stat<-"output_stat_TCGA"

library(future.apply)
library(future)

plan(multisession, workers = 25)

p_list <- list()
cutoff_list <- list()
coef_list <- list()

futures <- future_lapply(1:length(gene_list), function(j) {
  output <- Cheng_generateSurvInputfromTCGA(gene_list[j, 1], cancerType, dataDir)
  result <- Cheng_generateKMplot(output, outFile = paste0(path_out_pdf, cancerType, "_", gene_list[j, 1]))
  log_rank_p <- as.matrix(result$logRankP)
  cut_off <- as.matrix(result$EXPcut)
  coef <- as.matrix(result$coef)

  
  message(paste("Gene", j, "processed"))

  
  return(list(log_rank_p = log_rank_p, cut_off = cut_off, coef = coef))
})


for (j in 1:length(gene_list)) {

  result <- value(futures[[j]])

  p_list[[j]] <- result$log_rank_p
  cutoff_list[[j]] <- result$cut_off
  coef_list[[j]] <- result$coef
}

final_result <- data.frame(gene = gene_list[, 1], coef = unlist(coef_list),
                           cutoff = unlist(cutoff_list), p = unlist(p_list))

write.table(final_result, file = paste0(path_out_stat, cancerType, "_KM_stat_result_1.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

library(clusterProfiler)
library(dplyr)
library(tidyr)
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)

setwd("./Codes")


path_TCGA<-"./Example_data/"
setwd(path_TCGA)
load("CPTAC_PDAC_tumor_exp.Rdata")#load symbol_exp (your expression profiles, rownames is gene symbol)
background<-as.matrix(rownames(CPTAC_PDAC_tumor_exp))

path_raw<-"./Example_data/"
setwd(path_raw)
int_gene<-as.matrix(read.csv("common_pro-oncogenic_signatures.txt",header=F,sep="\t"))

pathway=enrichGO(gene=int_gene,OrgDb='org.Hs.eg.db',keyType = "SYMBOL",ont="BP",universe = background ,pAdjustMethod = "BH",qvalueCutoff=0.05)
path_table_all<-deg_GoTerm_clusterProfiler(pathway)
setwd(path_raw)
write.table(path_table_all,file="PDAC_GO_terms_KIRC_markers_pro-oncogenic_signatures.txt",sep="\t",row.names=F,col.names=T,quote=F)
