binary_metrics <- function(predictions, reference, numeric = TRUE, class_of_interest = 1)
{
  if(numeric)
  {
    if(class_of_interest == 1)
      other_class = 2
    else
      other_class = 1
    if(class_of_interest != 1 & class_of_interest !=2)
      print("incorrect parameter for class_of_interest")
    
    tp <- sum((predictions == class_of_interest) & (reference == class_of_interest))
    fp <- sum((predictions == class_of_interest) & (reference == other_class))
    fn <- sum((predictions == other_class) & (reference == class_of_interest))
    tn <- sum((predictions == other_class) & (reference == other_class))
    
    se <- tp/(tp+fn)
    sp <- tn/(tn+fp)
    ppv <- tp/(tp+fp)
    npv <- tn/(tn+fn)
    
    if(length(unique(predictions)) == 1)
    {
      area <- 0.50
    } else {
      if(class_of_interest == 1)
      {
        predictions_auc <- predictions
        predictions_auc[predictions_auc == 2] <- 0
        reference_auc <- reference
        reference_auc[reference_auc == 2] <- 0
        area <- auc(predictions_auc, reference)
      } else if(class_of_interest == 2)
      {
        predictions_auc <- predictions
        predictions_auc[predictions_auc == 1] <- 0
        reference_auc <- reference
        reference_auc[reference_auc == 1] <- 0
        predictions_auc[predictions_auc == 2] <- 1
        reference_auc[reference_auc == 2] <- 1
        area <- auc(predictions_auc, reference)
      }
    }
    metrics <- data.frame(se = round(se, 2), sp = round(sp, 2), ppv = round(ppv, 2), npv = round(npv, 2), auc = round(area, 2))
  }
  return (metrics)
}