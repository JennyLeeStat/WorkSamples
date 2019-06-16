library(tidyverse)
theme_set(theme_bw())

# read in the dataset
demo <- read.csv('data/roc_demo.csv')
comp <- read.csv('data/roc_comp.csv')
demo_comp <- read.csv('data/roc_demo_comp.csv')
pcs <- read.csv('data/roc_pcs.csv')
pcs_demo <- read.csv('data/roc_pcs_demo.csv')
pcs_demo_on_test <- read.csv('data/roc_best_on_test.csv')
densenet_on_test <- read.csv('data/roc_densenet.csv')



random_line <- data.frame(fpr = seq(0, 1, .01), 
                          tpr = seq(0, 1, .01),
                          name = 'random')

change_colnames <- function(df){
  df <- df[, 2:3]
  colnames(df) <- c('fpr', 'tpr')
  return(df)
}

demo <- change_colnames(demo)
comp <- change_colnames(comp)
demo_comp <- change_colnames(demo_comp)
pcs <- change_colnames(pcs)
pcs_demo <- change_colnames(pcs_demo)

# compare the performance of five shallow classifiers ======================
# Let's make one long data frame for the ROC figure

demo$Features <- 'demographic_only (AUC = 0.6637)'
comp$Features <- 'average T2 values (AUC = 0.5686)'
demo_comp$Features <- 'demographic + average T2 (AUC = 0.7047)'
pcs$Features <- 'PCs (AUC = 0.7532)'
pcs_demo$Features <- 'demographic + PCs (AUC = 0.7777)'

plot_df <- rbind(demo, comp, demo_comp, pcs, pcs_demo)
plot_df$Features <- factor(plot_df$Features, 
                       levels = c('average T2 values (AUC = 0.5686)',
                                  'demographic_only (AUC = 0.6637)',
                                  'demographic + average T2 (AUC = 0.7047)',
                                  'PCs (AUC = 0.7532)',
                                  'demographic + PCs (AUC = 0.7777)'))

ggplot(data = plot_df, 
       aes(x = fpr, y = tpr, linetype = Features)) +
  geom_line(colour = 'slategrey') + 
  # scale_colour_manual(
  #   values = c('slategrey', 'slateblue1', 'springgreen4',
  #              'tan2', 'tomato')) +
  xlab('1 - Specificity') +
  ylab('Sensitivity') + 
  theme(legend.position = c(0.7, 0.2),
        legend.background = element_blank()) +
  ggsave('figs/compare_aucs.png', width = 6, height = 6, dpi = 900)


# compare best two on the test set ==================================
densenet_on_test <- change_colnames(densenet_on_test)
pcs_demo_on_test <- change_colnames(pcs_demo_on_test)
densenet_on_test$Features <- 'DenseNet (AUC = 0.8244)'
pcs_demo_on_test$Features <- 'RF + demographic + PCs (AUC = 0.7777)'

plot_df2 <- rbind(densenet_on_test, pcs_demo_on_test)

ggplot(data = plot_df2, 
       aes(x = fpr, y = tpr, linetype = Features)) +
  geom_line(colour = 'slategrey') +
  theme(legend.position = c(0.7, 0.2)) +
  xlab('1 - Specificity') +
  ylab('Sensitivity') 


