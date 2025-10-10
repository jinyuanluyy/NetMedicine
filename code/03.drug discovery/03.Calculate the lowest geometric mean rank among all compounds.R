library(pheatmap)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(plyr)

cancer <- "PDAC"
load('16_landmark_new.Rdata')
load("04_CPTAC_CKAP2L_LINCS_predictions.RData")
load("02_landmark_land_drug_2020.Rdata")
landmark <- landmark_new
lincs_signatures <- t(landmark2020)
text <- lincs_signatures[1:100,1:100]
lincs_experiments <- land_drug_2020
table(lincs_experiments$cell_id)

lincs_experiments <- lincs_experiments %>%
  filter(pert_idose == "10 uM" & pert_itime == "24 h")

drug_preds <- results[[1]]
dz_sig <- results[[2]]

drug_preds$exp_id <- as.numeric(as.character(drug_preds$exp_id))
drug_preds_sig <- merge(drug_preds, lincs_experiments, by.x="exp_id", by.y="id")

CL <- "YAPC"

lincs_predictions_all <- subset(drug_preds_sig, cell_id %in% CL )

lincs_predictions_all <- lincs_predictions_all[order(lincs_predictions_all$cmap_score),]
table(lincs_predictions_all$cell_id)

control <- as.numeric(tail(lincs_predictions_all[order(lincs_predictions_all$cmap_score),], 1)$exp_id)


lincs_predictions <- lincs_predictions_all

drug_min_score <- aggregate(cmap_score ~ pert_iname, lincs_predictions, min)


drug_min_score$CPTAC_Rank <- rank(drug_min_score$cmap_score, ties.method = "min")
library(dplyr)


drug_min_score <- drug_min_score %>%
  left_join(lincs_predictions %>% dplyr::select(pert_iname, p), by = "pert_iname")

drug_min_score <- drug_min_score %>%
  distinct(pert_iname, .keep_all = TRUE)
CPTAC_CKAP2L_drug_score=drug_min_score
colnames(CPTAC_CKAP2L_drug_score)[4] <- "CPTAC_P"
colnames(CPTAC_CKAP2L_drug_score)[2] <- "CPTAC_cmap_score"
save(CPTAC_CKAP2L_drug_score,file = "CPTAC_CKAP2L_drug_score.Rdata")
