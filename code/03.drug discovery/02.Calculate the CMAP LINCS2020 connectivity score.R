subset_comparison_id <- "PDAC"
analysis_id <-  "cmap"
cmap_score <- function(sig_up, sig_down, drug_signature) {

  num_genes <- nrow(drug_signature)
  ks_up <- 0
  ks_down <- 0
  connectivity_score <- 0

  # I think we are re-ranking because the GeneID mapping changed the original rank range
  drug_signature[,"rank"] <- rank(drug_signature[,"rank"])

  # Merge the drug signature with the disease signature by GeneID. This becomes the V(j) from the algorithm description
  up_tags_rank <- merge(drug_signature, sig_up, by.x = "ids", by.y = 1)
  down_tags_rank <- merge(drug_signature, sig_down, by.x = "ids", by.y = 1)

  up_tags_position <- sort(up_tags_rank$rank)
  down_tags_position <- sort(down_tags_rank$rank)

  num_tags_up <- length(up_tags_position)
  num_tags_down <- length(down_tags_position)

  if(num_tags_up > 1 && num_tags_down > 1) {
    a_up <- 0
    b_up <- 0

    a_up <- max(sapply(1:num_tags_up,function(j) {
      j/num_tags_up - up_tags_position[j]/num_genes
    }))
    b_up <- max(sapply(1:num_tags_up,function(j) {
      up_tags_position[j]/num_genes - (j-1)/num_tags_up
    }))

    if(a_up > b_up) {
      ks_up <- a_up
    } else {
      ks_up <- -b_up
    }

    a_down <- 0
    b_down <- 0

    a_down <- max(sapply(1:num_tags_down,function(j) {
      j/num_tags_down - down_tags_position[j]/num_genes
    }))
    b_down <- max(sapply(1:num_tags_down,function(j) {
      down_tags_position[j]/num_genes - (j-1)/num_tags_down
    }))

    if(a_down > b_down) {
      ks_down <- a_down
    } else {
      ks_down <- -b_down
    }

    if (sum(sign(c(ks_down,ks_up))) == 0) {
      connectivity_score <- ks_up - ks_down # different signs
    }
  }
  return(connectivity_score)
}


landmark <- "1"
load("16_landmark_new.Rdata")
load("02_landmark_land_drug_2020.Rdata")
landmark_data <- landmark_new
gene.list <- landmark_data$gene_id
print(length(gene.list))
lincs_signatures <- t(landmark2020)
text <- lincs_signatures[1:100,1:100]

lincs_sig_info <- land_drug_2020
lincs_sig_info <- subset(lincs_sig_info, id %in% colnames(lincs_signatures))
#remove duplicate instances
lincs_sig_info <- lincs_sig_info[!duplicated(lincs_sig_info$id),]

sig.ids <- lincs_sig_info$id

# Load in the signature for the subset_comparison_id
dz_signature <- read.table("01_CPTAC_CKAP2L_cmap.txt",header=T,sep="\t")

dz_signature <- subset(dz_signature, GeneID %in% gene.list)

dz_genes_up <- subset(dz_signature,up_down=="up",select="GeneID")
dz_genes_down <- subset(dz_signature,up_down=="down",select="GeneID")
max_gene_size <- 250
#only select 100 genes
if (nrow(dz_genes_up)> max_gene_size){
  dz_genes_up <- data.frame(GeneID=dz_genes_up[1:max_gene_size,])
}
if (nrow(dz_genes_down)> max_gene_size){
  dz_genes_down <- data.frame(GeneID=dz_genes_down[1:max_gene_size,])
}

print ((nrow(dz_genes_up)+nrow(dz_genes_down)))


N_PERMUTATIONS <- 10000 #default 100000
random_sig_ids <- sample(colnames(lincs_signatures),N_PERMUTATIONS,replace=T)
count <- 0
random_cmap_scores <- NULL
for (exp_id in random_sig_ids){
  count <- count + 1
  print(count)
  if (landmark ==1){
    cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))
  }else{
    cmap_exp_signature <- data.frame(gene.list,  lincs_signatures[, as.character(exp_id)])
  }
  colnames(cmap_exp_signature) <- c("ids","rank")

  random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
  rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
  rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
  random_cmap_scores <- c(random_cmap_scores, cmap_score(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature))
}

rand_cmap_scores <- random_cmap_scores
save(rand_cmap_scores,file='03_CPTAC_CKAP2L_lincs_randoms.RData')


landmark <- "1"
load("16_landmark_new.Rdata")
load("02_landmark_land_drug_2020.Rdata")
landmark_data <- landmark_new
gene.list <- landmark_data$gene_id
lincs_signatures <- t(landmark2020)
text <- lincs_signatures[1:100,1:100]

lincs_sig_info <- land_drug_2020
lincs_sig_info <- subset(lincs_sig_info, id %in% colnames(lincs_signatures))
#remove duplicate instances
lincs_sig_info <- lincs_sig_info[!duplicated(lincs_sig_info$id),]

sig.ids <- lincs_sig_info$id
#duplicated_cells <- duplicated(lincs_sig_info$cell_id)

dz_signature <- read.table("01_CPTAC_CKAP2L_cmap.txt",header=T,sep="\t")

dz_signature <- subset(dz_signature, GeneID %in% gene.list)
dz_genes_up <- subset(dz_signature,up_down=="up",select="GeneID")
dz_genes_down <- subset(dz_signature,up_down=="down",select="GeneID")

max_gene_size <- 150
if (nrow(dz_genes_up)> max_gene_size){
  dz_genes_up <- data.frame(GeneID= dz_genes_up[1:max_gene_size,])
}
if (nrow(dz_genes_down)> max_gene_size){
  dz_genes_down <- data.frame(GeneID=dz_genes_down[1:max_gene_size,])
}

dz_cmap_scores <- NULL
count <- 0
for (exp_id in sig.ids) {
  count <- count + 1
  print(count)
  #print(paste("Computing score for disease against cmap_experiment_id =",exp_id))
  if (landmark ==1){
    cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))
  }else{
    cmap_exp_signature <- data.frame(gene.list,  lincs_signatures[, as.character(exp_id)])
  }
  colnames(cmap_exp_signature) <- c("ids","rank")
  dz_cmap_scores <- c(dz_cmap_scores, cmap_score(dz_genes_up,dz_genes_down,cmap_exp_signature))
}



load( '03_CPTAC_CKAP2L_cmap_randoms.RData')
random_scores <- unlist(rand_cmap_scores)
# Frequency-based p-value using absolute scores from sampling distribution to approximate two-tailed p-value
print("COMPUTING p-values")
p_values <- sapply(dz_cmap_scores,function(score) {
  length(which(abs(random_scores) >= abs(score))) / length(random_scores)
})
print("COMPUTING q-values")
q_values <- qvalue(p_values)$qvalues

drugs <- data.frame(exp_id = sig.ids, cmap_score = dz_cmap_scores, p = p_values, q = q_values, subset_comparison_id, analysis_id)

results = list(drugs, dz_signature)
save(results, file="04_CPTAC_CKAP2L_LINCS_predictions.RData")