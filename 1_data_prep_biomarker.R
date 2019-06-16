# =================================================================
# Filename: data_prep_biomarker.R
# Tasks: Detect and transform outliers for biomarker dataset
#     and build the dataset for the later analysis
# Author: Jenny Lee (JennyLee.Stat@gmail.com)
# Last updated: 1/24/2019
# ==================================================================


library(tidyverse)
library(readxl)

# clean up baseline biomarkers =====================================
 
# read in datasets 
bm_baseline <- read.csv(file = 'Data/biomarkers_baseline.csv')
bm_baseline$patient_id <- as.character(bm_baseline$patient_id)
bm_baseline <- bm_baseline[complete.cases(bm_baseline), ]


# Winsorization for outliers 
source('utils.R')

# Center the data before the principal component analysis

bm_scaled <- apply(bm_baseline[, -1], 2, function (x){
  # log transformation
  x = log(x)
  
  # winsorization
  x = winsorize(x)
  
  # standardization
  x = scale(x, center = TRUE, scale = FALSE)
  
  # # min-max scaling (between [0, 1])
  # x = (x - min(x)) / (max(x) - min(x))
  
  return(x)
})


bm_scaled <- data.frame(bm_scaled)
bm_scaled$patient_id <- as.character(bm_baseline$patient_id)
write.csv(bm_scaled, 'Data/biomarker_baseline_scaled.csv')


# Quick sanity check: Boxplot of transformed/scaled biomarker data
colnames(bm_scaled)
bm_scaled %>%
  select(-'patient_id') %>%
  boxplot(las=2, main = 'Log transformed and centered')


# build image dataset ==============================================
img_bl <- read.csv('Data/baseline_ipsi_image.csv')
img_bl$patient_id <- as.character(img_bl$patient_id)
img_bl <- img_bl[complete.cases(img_bl), ]

img_1yr <- read.csv('Data/oneyr_ipsi_image.csv')
img_1yr$patient_id <- as.character(img_1yr$patient_id)
img_1yr <- img_1yr[complete.cases(img_1yr), ]

ids = base::intersect(img_1yr$patient_id, img_bl$patient_id)
ids = base::intersect(ids, bm_scaled$patient_id)

img_bl <- img_bl %>% filter(patient_id %in% ids)
img_1yr <- img_1yr %>% filter(patient_id %in% ids)
bm_scaled <- bm_scaled %>% filter(patient_id %in% ids)


write.csv(ids, 'Data/patient_ids.csv')
write.csv(img_bl, 'Data/img_bl_small.csv')
write.csv(img_1yr, 'Data/img_1yr_small.csv')


# read in the clinical features dataset =============================
clinical <- read.csv('Data/clinical_features.csv')

clinical <- clinical[, c(2, 5, 15, 21, 22)]
clinical$patient_id <- as.character(clinical$patient_id)
clinical$effusion[is.na(clinical$effusion)] <- 0
clinical$cartilage_defects[is.na(clinical$cartilage_defects)] <- 0
clinical$synovitis_factor[is.na(clinical$synovitis_factor)] <- 0
clinical$meniscal_tear_severity[is.na(clinical$meniscal_tear_severity)] <- 0

# centered <- apply(clinical[, -1], 2, function(x){ 
#   x = scale(x, center = TRUE, scale = FALSE)
#   return(x)})
# 
# centered <- as.data.frame(centered)
# centered$patient_id <- as.character(clinical$patient_id)

combined <- bm_scaled %>% 
  left_join(clinical, by = 'patient_id') %>%
  filter(patient_id %in% ids)


write.csv(combined, file = 'Data/bm_scaled_clinical.csv')


# read in required datasets and merge them =========================

days_since_injury <- read.csv('Data/days_since_injury.csv')
days_since_injury <- days_since_injury[, c(1, 4)]
colnames(days_since_injury) <- c('patient_id', 'days_since_injury')
days_since_injury$patient_id <- as.character(days_since_injury$patient_id)


days_since_injury$days_since_injury <- scale(
  days_since_injury$days_since_injury, T, T)

demo <- read.csv('Data/demo_cleaned.csv')
demo$patient_id <- as.character(demo$patient_id)
demo <- demo[, c(1, 3, 6)]

demo[, -1] <- apply(demo[, -1], 2, function(x){
  scale(x, center = T, scale = T)})

combined <- left_join(combined, demo, by = 'patient_id')
combined <- left_join(combined, days_since_injury, by = 'patient_id')
combined <- combined %>% filter(patient_id %in% ids$patient_id)
combined$days_since_injury <- as.vector(combined$days_since_injury)

write.csv(combined, 'Data/combined_N25.csv')
