rm(list = ls()) 
library(ggplot2)
library(reshape2)

#set output path for generated figures, used only if output_to_pdf is TRUE
figuresOutputPath <- "/Users/Marco/Dropbox/"
output_to_pdf <- TRUE

#read data
tpehgdb_features <- read.csv("~/Dropbox/R workspace/github/physionet_oversampling/tpehgdb_features.csv")
#select third channel only, supposedly the best one for this analysis
tpehgdb_features <- tpehgdb_features[tpehgdb_features$channel == 3, ]

#Density plot for EHG features with respect to the outcome class (term/preterm)
selected_features <- c("rms", "fmed", "fpeak", "sample_entropy")
if(output_to_pdf)
{
  pdf(paste(figuresOutputPath,"figEDAden1.pdf", sep=""), width=8, height=7)#10 legend 8 no legend
}
p1 <- ggplot(data_to_plot, aes(value, fill = preterm)) + 
  geom_density(alpha = 0.6) +
  facet_wrap(~ variable, scales = "free") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  #xlab("condition")+
  ggtitle("Density plots for features extracted on the whole EHG recording")
p1
if(output_to_pdf)
{
  dev.off()
}

#leave one participant out cross-validation using all features, four classifiers and imbalanced data




