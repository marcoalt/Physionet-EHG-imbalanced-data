#clean workspace
rm(list = ls()) 
library(ggplot2) 
library(reshape2) 
library(tree)
library(e1071)
library(randomForest)
library(pROC)
library(xtable)
library(DMwR) #sampling

#set working directory to the current source file directory
dir.wd <- "Dropbox/R workspace/github/physionet_oversampling/"
source(paste(dir.wd, "/binary_metrics.R", sep = ""))
       
#set output path for generated figures, used only if output_to_pdf is TRUE
figuresOutputPath <- dir.wd
output_to_pdf <- TRUE

#read data
tpehgdb_features <- read.csv(paste(dir.wd, "/tpehgdb_features.csv", sep = ""))
#select third channel only, supposedly the best one for this analysis
tpehgdb_features <- tpehgdb_features[tpehgdb_features$channel == 3, ]

#Density plot for EHG features with respect to the outcome class (term/preterm)
selected_features <- c("rms", "fmed", "fpeak", "sample_entropy")
data_to_plot <- tpehgdb_features[, c("preterm", selected_features)]
data_to_plot <- melt(data_to_plot, id.vars = "preterm", measure.vars = selected_features)
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

#1: leave one participant out cross-validation using all features, four classifiers and imbalanced data
data_to_use <- tpehgdb_features
features <- c("rms", "fmed", "fpeak", "sample_entropy")

metrics_all <- data.frame()

#leave one participant out cross-validation
results_lr <- rep(NA, nrow(data_to_use))
results_tree <- rep(NA, nrow(data_to_use))
results_svm <- rep(NA, nrow(data_to_use))
results_rf <- rep(NA, nrow(data_to_use))

for(index_subj  in 1:nrow(data_to_use))
{
  #remove subject to validate
  training_data <- data_to_use[-index_subj, ]
  training_data_formula <- training_data[, c("preterm", features)]
  
  #select features in the validation set
  validation_data <- data_to_use[index_subj, features]
  
  #logistic regression
  glm.fit <- glm(preterm ~.,
                 data = training_data_formula,
                 family = binomial)
  glm.probs <- predict(glm.fit, validation_data, type = "response")
  predictions_lr <- ifelse(glm.probs < 0.5, "t", "f")
  results_lr[index_subj] <- predictions_lr
  
  #classification tree
  tree.fit <- tree(preterm ~.,
                   data = training_data_formula)
  predictions_tree <- predict(tree.fit, validation_data, type = "class")
  results_tree[index_subj] <- predictions_tree
  
  #svm
  svm <- svm(preterm ~.,
             data = training_data_formula
  )
  predictions_svm <- predict(svm, validation_data)
  results_svm[index_subj] <- predictions_svm
  
  #random forest      
  rf <- randomForest(preterm ~.,
                     data = training_data_formula)
  predictions_rf <- predict(rf, validation_data)
  results_rf[index_subj] <- predictions_rf   
}

#compute performance metrics
metrics_lr <- data.frame(binary_metrics(as.numeric(as.factor(results_lr)), as.numeric(data_to_use$preterm), class_of_interest = 2))
metrics_lr[, c("classifier")] <- c("logistic_regression")
metrics_all <- rbind(metrics_all, metrics_lr)

metrics_tree <- data.frame(binary_metrics(results_tree, as.numeric(data_to_use$preterm), class_of_interest = 2))
metrics_tree[, c("classifier")] <- c("tree")
metrics_all <- rbind(metrics_all, metrics_tree)

metrics_svm <- data.frame(binary_metrics(results_svm, as.numeric(data_to_use$preterm), class_of_interest = 2))
metrics_svm[, c("classifier")] <- c("svm")
metrics_all <- rbind(metrics_all, metrics_svm)

metrics_rf <- data.frame(binary_metrics(results_rf, as.numeric(data_to_use$preterm), class_of_interest = 2))
metrics_rf[, c("classifier")] <- c("random_forests")
metrics_all <- rbind(metrics_all, metrics_rf)  

if(output_to_pdf)
{
  pdf(paste(figuresOutputPath,"figResultsSeIgnore.pdf", sep=""), width=8, height=7)#10 legend 8 no legend
}
p1 <- ggplot(metrics_all, aes (x = classifier, y = se, fill = classifier)) + 
  geom_bar(stat = "identity") + 
  scale_y_continuous("Se", limits=c(0,1)) +
  scale_x_discrete("Classifier") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.line=element_blank(),axis.text.x=element_blank()) +
  ggtitle("Imbalanced data") +
  scale_fill_brewer(palette="Pastel2") +
  scale_color_brewer(palette="Pastel2")
p1
if(output_to_pdf)
{
  dev.off()
}
if(output_to_pdf)
{
  pdf(paste(figuresOutputPath,"figResultsSpIgnore.pdf", sep=""), width=8, height=7)#10 legend 8 no legend
}
p2 <- ggplot(metrics_all, aes (x = classifier, y = sp, fill = classifier)) + 
  geom_bar(stat = "identity") + 
  scale_y_continuous("Sp", limits=c(0,1)) +
  scale_x_discrete("Classifier") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.line=element_blank(),axis.text.x=element_blank()) +
  ggtitle("Imbalanced data") +
  scale_fill_brewer(palette="Pastel2") +
  scale_color_brewer(palette="Pastel2")
p2
if(output_to_pdf)
{
  dev.off()
}

#to_print <- metrics_all[, c("se", "sp", "auc", "classifier")]
#design.table <- xtable(to_print)
#print(design.table, floating=FALSE, include.rownames=FALSE )

#2: leave one participant out cross-validation using all features, four classifiers and undersampled data
data_to_use <- tpehgdb_features

metrics_all <- data.frame()

#leave one participant out cross-validation
results_lr <- rep(NA, nrow(data_to_use))
results_tree <- rep(NA, nrow(data_to_use))
results_svm <- rep(NA, nrow(data_to_use))
results_rf <- rep(NA, nrow(data_to_use))

rows_preterm <- sum(data_to_use$preterm == " t         ") #weird string, haven't changed it for now
for(index_subj  in 1:nrow(data_to_use))
{
  #remove subject to validate
  training_data <- data_to_use[-index_subj, ]
  training_data_preterm <- training_data[training_data$preterm == " t         ", ]
  training_data_term <- training_data[training_data$preterm == " f         ", ] 
  
  #get subsample to balance dataset
  indices <- sample(nrow(training_data_term), rows_preterm)
  training_data_term <- training_data_term[indices, ]
  training_data <- rbind(training_data_preterm, training_data_term)
  
  #select features in the training set
  training_data_formula <- training_data[, c("preterm", features)]
  
  #select features in the validation set
  validation_data <- data_to_use[index_subj, features]
  
  #logistic regression
  glm.fit <- glm(preterm ~.,
                 data = training_data_formula,
                 family = binomial)
  glm.probs <- predict(glm.fit, validation_data, type = "response")
  predictions_lr <- ifelse(glm.probs < 0.5, "t", "f")
  results_lr[index_subj] <- predictions_lr
  
  #classification tree
  tree.fit <- tree(preterm ~.,
                   data = training_data_formula)
  predictions_tree <- predict(tree.fit, validation_data, type = "class")
  results_tree[index_subj] <- predictions_tree
  
  #svm
  svm <- svm(preterm ~.,
                    data = training_data_formula
  )
  predictions_svm <- predict(svm, validation_data)
  results_svm[index_subj] <- predictions_svm

  #random forest      
  rf <- randomForest(preterm ~.,
                                     data = training_data_formula,
                                     sampsize = c(nrow(training_data_preterm), nrow(training_data_preterm)))
  predictions_rf <- predict(rf, validation_data)
  results_rf[index_subj] <- predictions_rf   
}

metrics_lr <- data.frame(binary_metrics(as.numeric(as.factor(results_lr)), as.numeric(data_to_use$preterm), class_of_interest = 2))
metrics_lr[, c("classifier")] <- c("logistic_regression")
metrics_all <- rbind(metrics_all, metrics_lr)

metrics_tree <- data.frame(binary_metrics(results_tree, as.numeric(data_to_use$preterm), class_of_interest = 2))
metrics_tree[, c("classifier")] <- c("tree")
metrics_all <- rbind(metrics_all, metrics_tree)

metrics_svm <- data.frame(binary_metrics(results_svm, as.numeric(data_to_use$preterm), class_of_interest = 2))
metrics_svm[, c("classifier")] <- c("svm")
metrics_all <- rbind(metrics_all, metrics_svm)

metrics_rf <- data.frame(binary_metrics(results_rf, as.numeric(data_to_use$preterm), class_of_interest = 2))
metrics_rf[, c("classifier")] <- c("random_forests")
metrics_all <- rbind(metrics_all, metrics_rf)  

if(output_to_pdf)
{
  pdf(paste(figuresOutputPath,"figResultsSeUndersampling.pdf", sep=""), width=8, height=7)#10 legend 8 no legend
}
p1 <- ggplot(metrics_all, aes (x = classifier, y = se, fill = classifier)) + 
  geom_bar(stat = "identity") + 
  scale_y_continuous("Se", limits=c(0,1)) +
  scale_x_discrete("Classifier") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.line=element_blank(),axis.text.x=element_blank()) +
  ggtitle("Undersampling") +
  scale_fill_brewer(palette="Pastel2") +
  scale_color_brewer(palette="Pastel2")
p1
if(output_to_pdf)
{
  dev.off()
}
if(output_to_pdf)
{
  pdf(paste(figuresOutputPath,"figResultsSpUndersampling.pdf", sep=""), width=8, height=7)#10 legend 8 no legend
}
p2 <- ggplot(metrics_all, aes (x = classifier, y = sp, fill = classifier)) + 
  geom_bar(stat = "identity") + 
  scale_y_continuous("Sp", limits=c(0,1)) +
  scale_x_discrete("Classifier") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.line=element_blank(),axis.text.x=element_blank()) +
  ggtitle("Undersampling") +
  scale_fill_brewer(palette="Pastel2") +
  scale_color_brewer(palette="Pastel2")
p2
if(output_to_pdf)
{
  dev.off()
}

#to_print <- metrics_all[, c("se", "sp", "auc", "classifier")]
#design.table <- xtable(to_print)
#print(design.table, floating=FALSE, include.rownames=FALSE )

#3: leave one participant out cross-validation using all features, four classifiers and oversampled data
#SMOTE is performed before cross-validation, i.e. this is bad cross-validation
data_to_use <- tpehgdb_features
data_to_use_smote <- SMOTE(preterm ~ . , cbind(data_to_use[, c("preterm", features)]), k=5, perc.over = 600)

metrics_all <- data.frame()

#leave one participant out cross-validation
results_lr <- rep(NA, nrow(data_to_use_smote))
results_tree <- rep(NA, nrow(data_to_use_smote))
results_svm <- rep(NA, nrow(data_to_use_smote))
results_rf <- rep(NA, nrow(data_to_use_smote))

for(index_subj  in 1:nrow(data_to_use_smote))
{
  #remove subject to validate
  training_data <- data_to_use_smote[-index_subj, ]
  
  #no need to balance the dataset anymore     
  #select features in the training set
  training_data_formula <- training_data[, c("preterm", features)]

  #select features in the validation set
  validation_data <- data_to_use_smote[index_subj, features]
  
  #logistic regression
  glm.fit <- glm(preterm ~.,
                 data = training_data_formula,
                 family = binomial)
  glm.probs <- predict(glm.fit, validation_data, type = "response")
  predictions_lr <- ifelse(glm.probs < 0.5, "t", "f")
  results_lr[index_subj] <- predictions_lr
  
  #classification tree
  tree.fit <- tree(preterm ~.,
                   data = training_data_formula)
  predictions_tree <- predict(tree.fit, validation_data, type = "class")
  results_tree[index_subj] <- predictions_tree
  
  #svm
  svm <- svm(preterm ~.,
             data = training_data_formula
  )
  predictions_svm <- predict(svm, validation_data)
  results_svm[index_subj] <- predictions_svm
  
  #random forest      
  rf <- randomForest(preterm ~.,
                     data = training_data_formula)
  predictions_rf <- predict(rf, validation_data)
  results_rf[index_subj] <- predictions_rf   
}

metrics_lr <- data.frame(binary_metrics(as.numeric(as.factor(results_lr)), as.numeric(data_to_use_smote$preterm), class_of_interest = 2))
metrics_lr[, c("classifier")] <- c("logistic_regression")
metrics_all <- rbind(metrics_all, metrics_lr)

metrics_tree <- data.frame(binary_metrics(results_tree, as.numeric(data_to_use_smote$preterm), class_of_interest = 2))
metrics_tree[, c("classifier")] <- c("tree")
metrics_all <- rbind(metrics_all, metrics_tree)

metrics_svm <- data.frame(binary_metrics(results_svm, as.numeric(data_to_use_smote$preterm), class_of_interest = 2))
metrics_svm[, c("classifier")] <- c("svm")
metrics_all <- rbind(metrics_all, metrics_svm)

metrics_rf <- data.frame(binary_metrics(results_rf, as.numeric(data_to_use_smote$preterm), class_of_interest = 2))
metrics_rf[, c("classifier")] <- c("random_forests")
metrics_all <- rbind(metrics_all, metrics_rf)  

if(output_to_pdf)
{
  pdf(paste(figuresOutputPath,"figResultsSeOversamplingBad.pdf", sep=""), width=8, height=7)#10 legend 8 no legend
}
p1 <- ggplot(metrics_all, aes (x = classifier, y = se, fill = classifier)) + 
  geom_bar(stat = "identity") + 
  scale_y_continuous("Se", limits=c(0,1)) +
  scale_x_discrete("Classifier") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.line=element_blank(),axis.text.x=element_blank()) +
  ggtitle("Oversampling done wrong") +
  scale_fill_brewer(palette="Pastel2") +
  scale_color_brewer(palette="Pastel2")
p1
if(output_to_pdf)
{
  dev.off()
}
if(output_to_pdf)
{
  pdf(paste(figuresOutputPath,"figResultsSpOversamplingBad.pdf", sep=""), width=8, height=7)#10 legend 8 no legend
}
p2 <- ggplot(metrics_all, aes (x = classifier, y = sp, fill = classifier)) + 
  geom_bar(stat = "identity") + 
  scale_y_continuous("Sp", limits=c(0,1)) +
  scale_x_discrete("Classifier") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.line=element_blank(),axis.text.x=element_blank()) +
  ggtitle("Oversampling done wrong") +
  scale_fill_brewer(palette="Pastel2") +
  scale_color_brewer(palette="Pastel2")
p2
if(output_to_pdf)
{
  dev.off()
}

#to_print <- metrics_all[, c("se", "sp", "auc", "classifier")]
#design.table <- xtable(to_print)
#print(design.table, floating=FALSE, include.rownames=FALSE )

#4: leave one participant out cross-validation using all features, four classifiers and oversampled data
#SMOTE is performed during cross-validation, i.e. this is good cross-validation
data_to_use <- tpehgdb_features

metrics_all <- data.frame()

#leave one participant out cross-validation
results_lr <- rep(NA, nrow(data_to_use))
results_tree <- rep(NA, nrow(data_to_use))
results_svm <- rep(NA, nrow(data_to_use))
results_rf <- rep(NA, nrow(data_to_use))

for(index_subj  in 1:nrow(data_to_use))
{
  #remove subject to validate
  training_data <- data_to_use[-index_subj, ]
  training_data_smote <- SMOTE(preterm ~ . , cbind(training_data[, c("preterm", features)]), k=5, perc.over = 600)
  
  #no need to balance the dataset anymore     
  #select features in the training set
  training_data_formula <- training_data_smote[, c("preterm", features)]
  
  #select features in the validation set
  validation_data <- data_to_use[index_subj, features]
  
  #logistic regression
  glm.fit <- glm(preterm ~.,
                 data = training_data_formula,
                 family = binomial)
  glm.probs <- predict(glm.fit, validation_data, type = "response")
  predictions_lr <- ifelse(glm.probs < 0.5, "t", "f")
  results_lr[index_subj] <- predictions_lr
  
  #classification tree
  tree.fit <- tree(preterm ~.,
                   data = training_data_formula)
  predictions_tree <- predict(tree.fit, validation_data, type = "class")
  results_tree[index_subj] <- predictions_tree
  
  #svm
  svm <- svm(preterm ~.,
             data = training_data_formula
  )
  predictions_svm <- predict(svm, validation_data)
  results_svm[index_subj] <- predictions_svm
  
  #random forest      
  rf <- randomForest(preterm ~.,
                     data = training_data_formula)
  predictions_rf <- predict(rf, validation_data)
  results_rf[index_subj] <- predictions_rf   
}

metrics_lr <- data.frame(binary_metrics(as.numeric(as.factor(results_lr)), as.numeric(data_to_use$preterm), class_of_interest = 2))
metrics_lr[, c("classifier")] <- c("logistic_regression")
metrics_all <- rbind(metrics_all, metrics_lr)

metrics_tree <- data.frame(binary_metrics(results_tree, as.numeric(data_to_use$preterm), class_of_interest = 2))
metrics_tree[, c("classifier")] <- c("tree")
metrics_all <- rbind(metrics_all, metrics_tree)

metrics_svm <- data.frame(binary_metrics(results_svm, as.numeric(data_to_use$preterm), class_of_interest = 2))
metrics_svm[, c("classifier")] <- c("svm")
metrics_all <- rbind(metrics_all, metrics_svm)

metrics_rf <- data.frame(binary_metrics(results_rf, as.numeric(data_to_use$preterm), class_of_interest = 2))
metrics_rf[, c("classifier")] <- c("random_forests")
metrics_all <- rbind(metrics_all, metrics_rf)  

if(output_to_pdf)
{
  pdf(paste(figuresOutputPath,"figResultsSeOversamplingGood.pdf", sep=""), width=8, height=7)#10 legend 8 no legend
}
p1 <- ggplot(metrics_all, aes (x = classifier, y = se, fill = classifier)) + 
  geom_bar(stat = "identity") + 
  scale_y_continuous("Se", limits=c(0,1)) +
  scale_x_discrete("Classifier") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.line=element_blank(),axis.text.x=element_blank()) +
  ggtitle("Oversampling done properly") +
  scale_fill_brewer(palette="Pastel2") +
  scale_color_brewer(palette="Pastel2")
p1
if(output_to_pdf)
{
  dev.off()
}
if(output_to_pdf)
{
  pdf(paste(figuresOutputPath,"figResultsSpOversamplingGood.pdf", sep=""), width=8, height=7)#10 legend 8 no legend
}
p2 <- ggplot(metrics_all, aes (x = classifier, y = sp, fill = classifier)) + 
  geom_bar(stat = "identity") + 
  scale_y_continuous("Sp", limits=c(0,1)) +
  scale_x_discrete("Classifier") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.line=element_blank(),axis.text.x=element_blank()) +
  ggtitle("Oversampling done properly") +
  scale_fill_brewer(palette="Pastel2") +
  scale_color_brewer(palette="Pastel2")
p2
if(output_to_pdf)
{
  dev.off()
}

#to_print <- metrics_all[, c("se", "sp", "auc", "classifier")]
#design.table <- xtable(to_print)
#print(design.table, floating=FALSE, include.rownames=FALSE )
