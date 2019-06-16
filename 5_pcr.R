# =================================================================
# Filename: 5_pcr.R
# Tasks: Preliminary analysis for data dimension reduction 
# Author: Jenny Lee (JennyLee.Stat@gmail.com)
# Last updated: 1/29/2019
# ==================================================================

library(tidyverse)
library(xtable)
library(vtreat)
library(lmtest)
library(pls)
library(MASS)
source('utils.R')
theme_set(theme_bw())

# read in the datasets =============================================
bm <- read.csv('Data/combined_N25.csv')
rownames(bm) <-  as.character(bm$patient_id)

# simplify meniscus tear feature to two level factor variable
bm$meniscal_tear_severity[bm$meniscal_tear_severity >= 1] <- 1

# let's make a smaller dataset with only continuous dataset
bm_all <- bm[, c(-1, -20)]
bm <- bm[, 2:19]

categorical <- c("effusion", "cartilage_defects",     
                 "synovitis_factor", "meniscal_tear_severity")
categorical <- bm_all[, categorical]
count(categorical, effusion, cartilage_defects,
      synovitis_factor, meniscal_tear_severity)


vbr_vol <- read.csv('Data/VBR_cluster_vols_by_compartments.csv')
rownames(vbr_vol) <- vbr_vol$patient_id
vol_mf <- vbr_vol[rownames(bm),  "vol_MF"]
vol_lf <- vbr_vol[rownames(bm),  "vol_LF"]
vol_mt <- vbr_vol[rownames(bm),  "vol_MT"]
vol_lt <- vbr_vol[rownames(bm),  "vol_LT"]
vol_tro <- vbr_vol[rownames(bm),  "vol_TRO"]
vol_pat <- vbr_vol[rownames(bm),  "vol_PAT"]

bm$vol_mf <- vol_mf
bm_all$vol_mf <- vol_mf

img_1yr <- read.csv('Data/img_1yr_small.csv')
t1rho_mf <- img_1yr$t1rho__global_mf

# y-aware scale =================================================
examplePruneSig <- 1.0
treatmentsM <- designTreatmentsN(
  bm_all,
  setdiff(colnames(bm_all), 'vol_mf'),
  'vol_mf',
  verbose = TRUE)

cont_data_scaled <- prepare(treatmentsM,
                            bm_all,
                            pruneSig = examplePruneSig,
                            # pruning is off by setting the sig as zero
                            scale = TRUE)

cont_data_scaled$vol_mf <- vol_mf
scaled_bm <- cont_data_scaled[, c( -26 )]
rownames(scaled_bm) <- rownames(bm_all)
boxplot(scaled_bm, las= 3)

# cont_data_scaled <- apply(scaled_bm, 2, function(x){
#   y_aware_scale(x, y = vol_mf)
# })
# 
# cont_data_scaled <- as.data.frame(cont_data_scaled)
# rownames(cont_data_scaled) <- rownames(bm)
# boxplot(cont_data_scaled, las = 2) # doesn't look like scaled at all

# run the PCR ===================================================

pca_res <- prcomp(scaled_bm, center = TRUE, scale = FALSE)
par(mfrow=c(1, 3))
biplot(pca_res, choices=1:2)
biplot(pca_res, choices=c(1, 3))
biplot(pca_res, choices=c(1, 4))


# plot the magnitude of singular values
tmp <- lapply(1:length(pca_res$sdev), function(x){
  paste('PC', x, sep='')
})

tmp_df <- data.frame(pc_name = unlist(tmp),
                     pc = 1:length(pca_res$sdev),
                     magnitude = (pca_res$sdev) ** 2,
                     stringsAsFactors = FALSE)


ggplot(tmp_df[1:12, ], 
       aes(x = pc, y = magnitude, ymax = magnitude)) +
  geom_point(size=2, colour="#993399") + 
  geom_linerange(aes(ymin = 0)) + 
  scale_x_continuous(breaks = 1:12, labels = tmp_df$pc_name[1:12]) + 
  xlab(' ') + 
  ylab('Variances') + 
  ggtitle(' ') +
  ggsave('Figs/scree_plot_scaled_bm.png', 
         width = 7, height = 3.5)

  
# plot the scaled loadings 
rot3 <- pca_res$rotation[, 1:4]
rot3 <- as.data.frame(rot3)
rot3$varName = rownames(rot3)

rot3_long <- rot3 %>% 
  gather(key = "PC", "loadings", -varName) 

ggplot(rot3_long, aes(x = varName, y = loadings)) +
  geom_point(color ="red", size = 5, shape = 18, alpha = .75) + 
  geom_linerange(aes(ymax = loadings, ymin = 0), 
                 alpha = .5, size = 2) +
  facet_wrap(~PC, nrow = 1) + 
  coord_flip() + 
  ylab(" ") +
  xlab(" ") +
  ggtitle("Scaled variable loadings (First 4 PCs)") +
  ggsave('Figs/loadings_scaled_bm.png', width = 8, height = 6)


# preliminary analysis 2: T1rho measurement in yr 1 as outcome ===========
bm_all <- bm_all[, -26]
bm_all$t1rho_mf <- t1rho_mf 
bm_scaled_global_mf <- get_scaled_df(bm_all, "t1rho_mf")

pca_res <- prcomp(bm_scaled_global_mf, center = TRUE, scale = FALSE)
projection <- pca_res$rotation[, 1:10]
projected_data <- as.data.frame(
  as.matrix(bm_scaled_global_mf) %*% projection,
  stringsAsFactors = FALSE)  
projected_data$t1rho_mf <- t1rho_mf
full_model <- lm(t1rho_mf ~ ., data = projected_data)
summary(full_model)
step_model <- stepAIC(full_model, direction = "backward", trace = FALSE)
step_model$anova
summary(step_model)



# Primary analysis: Vol MF as outcome =====================================

# pcr_df <- scaled_bm
# pcr_df$vol_mf <- vol_mf
# mod <- pcr(vol_mf ~. , ncomp = 20, data = pcr_df, scale = FALSE)
# plot(mod, labels = rownames(scaled_bm))
# summary(mod)

pca_res <- prcomp(scaled_bm, center = TRUE, scale = FALSE)
projection <- pca_res$rotation[, 1:10]
projected_data <- as.data.frame(
  as.matrix(scaled_bm) %*% projection,
  stringsAsFactors = FALSE)  

projected_data$vol_mf <- vol_mf

# stepwise selection
full_model <- lm(vol_mf ~ ., data = projected_data)
step_model <- stepAIC(full_model, direction = "backward", trace = FALSE)
step_model$anova

summary(full_model)
summary(step_model)

res = summary(step_model)$coefficients
V = pca_res$rotation[, c(1, 2, 10)]
beta = as.matrix(res[2:4, 1], nrow=3)
beta_lower = as.matrix(res[2:4, 1] - res[2:4, 2], nrow=3)
beta_upper = as.matrix(res[2:4, 1] + res[2:4, 2], nrow=3)
new_beta = V %*% beta
new_beta_lower = V %*% beta_lower
new_beta_upper = V %*% beta_upper
recovered = as.data.frame(cbind(new_beta, new_beta_lower, new_beta_upper))
colnames(recovered) = c('Estimate', 'CI_lower', 'CI_upper')
tmp_l = pmin(recovered$CI_lower, recovered$CI_upper)
tmp_u = pmax(recovered$CI_lower, recovered$CI_upper)
recovered$CI_lower <- tmp_l
recovered$CI_upper <- tmp_u
colnames(recovered) = c('Estimate', 'CI_lower', 'CI_upper')
recovered



# weights plot for significant variates

rot <- pca_res$rotation[, c(2, 3, 5, 7)] %>%
  as.data.frame() %>%
  mutate(varName = rownames(pca_res$rotation)) 
rownames(rot) <- rot[, 5]
rot[, 1:5] <- round(rot[, 1:4], 4)
rot <- rot[, -5]
rot
xtable(rot)

# plot the scaled loadings 

rot3 <- as.data.frame(rot)
rot3$varName = rownames(rot3)


rot3_long <- rot3 %>% 
  gather(key = "PC", "loadings", -varName) 


ggplot(rot3_long, aes(x = varName, y = loadings)) +
  geom_point(color = 'seagreen4', size = 3) + 
  geom_linerange(aes(ymax = loadings, ymin = 0), 
                 color = 'seagreen4') +
  facet_wrap(~PC, nrow = 1) + 
  coord_flip() + 
  ylab(" ") +
  xlab(" ") + 
  ggtitle(" ") + 
  ggsave("Figs/sig_coef.png", width = 8, height = 6)


# exploratory analysis: vol LF as outcome =============================

bm_all <- bm_all[, -26]
bm_all$vol_lf <- vol_lf
bm_scaled_lf <- get_scaled_df(bm_all, "vol_lf")

pca_res_lf <- prcomp(bm_scaled_lf, center = TRUE, scale = FALSE)
projection <- pca_res_lf$rotation[, 1:10]
projected_data <- as.data.frame(
  as.matrix(bm_scaled_lf) %*% projection,
  stringsAsFactors = FALSE)  

projected_data$vol_lf <- vol_lf

# stepwise selection
full_model <- lm(vol_lf ~ ., data = projected_data)
step_model <- stepAIC(full_model, direction = "backward", trace = FALSE)
step_model$anova

summary(full_model)
summary(step_model)

res = summary(step_model)$coefficients
V = pca_res_lf$rotation[, c(2, 5, 8, 10)]
beta = as.matrix(res[2:5, 1], nrow=3)
beta_lower = as.matrix(res[2:5, 1] - res[2:5, 2], nrow=3)
beta_upper = as.matrix(res[2:5, 1] + res[2:5, 2], nrow=3)
new_beta = V %*% beta
new_beta_lower = V %*% beta_lower
new_beta_upper = V %*% beta_upper
recovered = as.data.frame(cbind(new_beta, new_beta_lower, new_beta_upper))
colnames(recovered) = c('Estimate', 'CI_lower', 'CI_upper')
tmp_l = pmin(recovered$CI_lower, recovered$CI_upper)
tmp_u = pmax(recovered$CI_lower, recovered$CI_upper)
recovered$CI_lower <- tmp_l
recovered$CI_upper <- tmp_u
colnames(recovered) = c('Estimate', 'CI_lower', 'CI_upper')
recovered

# exploratory analysis: vol LT as outcome ===============================
bm_all <- bm_all[, -26]
bm_all$vol_lt <- vol_lt
bm_scaled_lt <- get_scaled_df(bm_all, "vol_lt")

pca_res_lt <- prcomp(bm_scaled_lt, center = TRUE, scale = FALSE)
projection <- pca_res_lt$rotation[, 1:10]
projected_data <- as.data.frame(
  as.matrix(bm_scaled_lt) %*% projection,
  stringsAsFactors = FALSE)  

projected_data$vol_lt <- vol_lt

# stepwise selection
full_model <- lm(vol_lt ~ ., data = projected_data)
step_model <- stepAIC(full_model, direction = "backward", trace = FALSE)
step_model$anova

summary(full_model)
summary(step_model)


res = summary(step_model)$coefficients
V = pca_res_lt$rotation[, c(1)]
beta = as.matrix(res[2, 1], nrow=3)
beta_lower = as.matrix(res[2, 1] - res[2, 2], nrow=3)
beta_upper = as.matrix(res[2, 1] + res[2, 2], nrow=3)
new_beta = V %*% beta
new_beta_lower = V %*% beta_lower
new_beta_upper = V %*% beta_upper
recovered = as.data.frame(cbind(new_beta, new_beta_lower, new_beta_upper))
colnames(recovered) = c('Estimate', 'CI_lower', 'CI_upper')
rownames(recovered) = names(V)
tmp_l = pmin(recovered$CI_lower, recovered$CI_upper)
tmp_u = pmax(recovered$CI_lower, recovered$CI_upper)
recovered$CI_lower <- tmp_l
recovered$CI_upper <- tmp_u
colnames(recovered) = c('Estimate', 'CI_lower', 'CI_upper')
recovered


# exploratory analysis: vol MT as outcome ===============================

bm_all <- bm_all[, -26]
bm_all$vol_mt <- vol_mt
bm_scaled_mt <- get_scaled_df(bm_all, "vol_mt")

pca_res_mt <- prcomp(bm_scaled_mt, center = TRUE, scale = FALSE)
projection <- pca_res_mt$rotation[, 1:10]
projected_data <- as.data.frame(
  as.matrix(bm_scaled_mt) %*% projection,
  stringsAsFactors = FALSE)  

projected_data$vol_mt <- vol_mt

# stepwise selection
full_model <- lm(vol_mt ~ ., data = projected_data)
step_model <- stepAIC(full_model, direction = "backward", trace = FALSE)
step_model$anova

summary(full_model)
summary(step_model)


res = summary(step_model)$coefficients
V = pca_res_mt$rotation[, c(1, 2, 3, 8, 9)]
beta = as.matrix(res[2:6, 1], nrow=3)
beta_lower = as.matrix(res[2:6, 1] - res[2:6, 2], nrow=3)
beta_upper = as.matrix(res[2:6, 1] + res[2:6, 2], nrow=3)
new_beta = V %*% beta
new_beta_lower = V %*% beta_lower
new_beta_upper = V %*% beta_upper
recovered = as.data.frame(cbind(new_beta, new_beta_lower, new_beta_upper))
colnames(recovered) = c('Estimate', 'CI_lower', 'CI_upper')
tmp_l = pmin(recovered$CI_lower, recovered$CI_upper)
tmp_u = pmax(recovered$CI_lower, recovered$CI_upper)
recovered$CI_lower <- tmp_l
recovered$CI_upper <- tmp_u
colnames(recovered) = c('Estimate', 'CI_lower', 'CI_upper')
recovered

# exploratory analysis: vol PAT as outcome ===============================

bm_all <- bm_all[, -26]
bm_all$vol_pat <- vol_pat
bm_scaled_pat <- get_scaled_df(bm_all, "vol_pat")

pca_res_pat <- prcomp(bm_scaled_pat, center = TRUE, scale = FALSE)
projection <- pca_res_pat$rotation[, 1:10]
projected_data <- as.data.frame(
  as.matrix(bm_scaled_pat) %*% projection,
  stringsAsFactors = FALSE)  

projected_data$vol_pat <- vol_pat

# stepwise selection
full_model <- lm(vol_pat ~ ., data = projected_data)
step_model <- stepAIC(full_model, direction = "backward", trace = FALSE)
step_model$anova

summary(full_model)
summary(step_model)


res = summary(step_model)$coefficients
V = pca_res_pat$rotation[, c(1, 2, 6)]
beta = as.matrix(res[2:4, 1], nrow=3)
beta_lower = as.matrix(res[2:4, 1] - res[2:4, 2], nrow=3)
beta_upper = as.matrix(res[2:4, 1] + res[2:4, 2], nrow=3)
new_beta = V %*% beta
new_beta_lower = V %*% beta_lower
new_beta_upper = V %*% beta_upper
recovered = as.data.frame(cbind(new_beta, new_beta_lower, new_beta_upper))
colnames(recovered) = c('Estimate', 'CI_lower', 'CI_upper')
tmp_l = pmin(recovered$CI_lower, recovered$CI_upper)
tmp_u = pmax(recovered$CI_lower, recovered$CI_upper)
recovered$CI_lower <- tmp_l
recovered$CI_upper <- tmp_u
colnames(recovered) = c('Estimate', 'CI_lower', 'CI_upper')
recovered

# exploratory analysis: vol TRO as outcome ===============================

bm_all <- bm_all[, -26]
bm_all$vol_tro <- vol_tro
bm_scaled_tro <- get_scaled_df(bm_all, "vol_tro")

pca_res_tro <- prcomp(bm_scaled_tro, center = TRUE, scale = FALSE)
projection <- pca_res_tro$rotation[, 1:10]
projected_data <- as.data.frame(
  as.matrix(bm_scaled_tro) %*% projection,
  stringsAsFactors = FALSE)  

projected_data$vol_tro <- vol_tro

# stepwise selection
full_model <- lm(vol_tro ~ ., data = projected_data)
step_model <- stepAIC(full_model, direction = "backward", trace = FALSE)
step_model$anova

summary(full_model)
summary(step_model)


res = summary(step_model)$coefficients
V = pca_res_pat$rotation[, c(2, 7)]
beta = as.matrix(res[2:3, 1], nrow=3)
beta_lower = as.matrix(res[2:3, 1] - res[2:3, 2], nrow=3)
beta_upper = as.matrix(res[2:3, 1] + res[2:3, 2], nrow=3)
new_beta = V %*% beta
new_beta_lower = V %*% beta_lower
new_beta_upper = V %*% beta_upper
recovered = as.data.frame(cbind(new_beta, new_beta_lower, new_beta_upper))
colnames(recovered) = c('Estimate', 'CI_lower', 'CI_upper')
tmp_l = pmin(recovered$CI_lower, recovered$CI_upper)
tmp_u = pmax(recovered$CI_lower, recovered$CI_upper)
recovered$CI_lower <- tmp_l
recovered$CI_upper <- tmp_u
colnames(recovered) = c('Estimate', 'CI_lower', 'CI_upper')
recovered
