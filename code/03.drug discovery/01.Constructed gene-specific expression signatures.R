library(limma)
library(dplyr)
library(tibble)


load("CPTAC_PDAC_tumor_exp.Rdata")

genes <- c("CKAP2L")
gene_data <- CPTAC_PDAC_tumor_exp[genes, ]

get_top_bottom_samples <- function(expression_values, top_percent = 0.1, bottom_percent = 0.1) {
  sorted_values <- sort(expression_values, decreasing = TRUE)
  n <- length(expression_values)
  top_n <- ceiling(n * top_percent)
  bottom_n <- ceiling(n * bottom_percent)
  top_samples <- names(sorted_values[1:top_n])
  bottom_samples <- names(sorted_values[(n - bottom_n + 1):n])
  list(top_samples = top_samples, bottom_samples = bottom_samples)
}

results <- lapply(genes, function(gene) {
  expression_values <- as.numeric(gene_data[gene, ])
  names(expression_values) <- colnames(gene_data)
  get_top_bottom_samples(expression_values)
})


extract_expression_profile <- function(samples) {
  CPTAC_PDAC_tumor_exp[, samples]
}


top_samples <- unique(unlist(lapply(results, function(result) result$top_samples)))
bottom_samples <- unique(unlist(lapply(results, function(result) result$bottom_samples)))


top_CPTAC_CKAP2L <- extract_expression_profile(top_samples)


bottom_CPTAC_CKAP2L <- extract_expression_profile(bottom_samples)

save(top_CPTAC_CKAP2L,bottom_CPTAC_CKAP2L,file = "./CPTAC_CKAP2L.RData")


load("CPTAC_PDAC_tumor_exp.Rdata")
load("CPTAC_CKAP2L.RData")

result1 <- bottom_CPTAC_CKAP2L
result2 <- top_CPTAC_CKAP2L
merged_df <- cbind(result1, result2)
group <- rep(c("con", "treat"), times = c(ncol(result1), ncol(result2)))
group <- factor(group, levels = c("con", "treat"))
design <- model.matrix(~group)
fit <- lmFit(merged_df, design)
fit <- eBayes(fit)
ALLDIFF <- topTable(fit, coef=2, number = Inf)


logFC_t <- 1
p_t <- 0.05

k1 <- (ALLDIFF$P.Value < p_t) & (ALLDIFF$logFC < -logFC_t)
k2 <- (ALLDIFF$P.Value < p_t) & (ALLDIFF$logFC > logFC_t)
diffgene1 <- mutate(ALLDIFF, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))


CPTAC_CKAP2L_logFC <- filter(diffgene1, change %in% c("up", "down"))
table(CPTAC_CKAP2L_logFC$change)

save(CPTAC_CKAP2L_logFC,file = "CPTAC_CKAP2L_logFC.Rdata")

cancer =  'PDAC' 

CPTAC_CKAP2L_logFC$id <- rownames(CPTAC_CKAP2L_logFC)
rownames(CPTAC_CKAP2L_logFC) <- NULL
dz_signature=CPTAC_CKAP2L_logFC
dz_signature <- subset(dz_signature, select=c("id", "logFC", "P.Value"))
names(dz_signature) <- c("GeneID",  "value", "pval")
dz_signature <- subset(dz_signature, !is.na(value))
dz_signature <- dz_signature[order(dz_signature$value),]

library("clusterProfiler")
library(org.Hs.eg.db)
library(dplyr)
library(stringr)


data = dz_signature$GeneID
entrez_id <- bitr(data,
                  fromType="SYMBOL", 
                  toType="ENTREZID", 
                  OrgDb="org.Hs.eg.db")
identical(entrez_id $SYMBOL,dz_signature$GeneID)
rownames(dz_signature) <- dz_signature$GeneID
dz_signature <- dz_signature[entrez_id $SYMBOL,]
dz_signature$GeneID <- entrez_id$ENTREZID
dz_signature <- dz_signature[order(dz_signature$value),]

dz_signature$up_down <- "up"
dz_signature$up_down[dz_signature$value < 0] <- "down"

write.table(dz_signature, "01_CPTAC_CKAP2L_cmap.txt", quote=F, col.names=T, row.names=F, sep="\t" )
