################################################################################
# Calculate optimal number of components from downsampling plsda model objects
################################################################################
getNComps = function(directory, ncomp) {
  require(matrixStats)
  models = list.files(path = directory)
  accuracies = matrix(nrow = ncomp)
  
  for(i in 1:length(models)) {
    model = readRDS(paste(directory, models[i], sep = "/"))
    accuracies = cbind(accuracies, as.matrix(model$results$Accuracy))
  }
  
  accuracies = accuracies[,-1]
  meanAcc = as.matrix(rowMeans(accuracies))
  sdAcc = as.matrix(rowSds(accuracies))
  lowerSd = meanAcc - sdAcc
  #get accuracies less than the standard deviation of the highest mean accuracy
  lowerMeans = meanAcc[meanAcc < lowerSd[which.max(meanAcc)]]
  optimalComp = length(lowerMeans) + 1 #first component with an average 
  # accuracy within 1 sd of the component with the highest mean accuracy
  
  return(optimalComp)
}

saveRDS(getNComps, "revised pls/functions/getNComps.rds")
