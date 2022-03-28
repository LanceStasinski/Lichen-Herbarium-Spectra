################################################################################
# PLS function
################################################################################

runPlsda = function(iteration, spectra, className, ncomp, resampling, include_age,
                    modelDirectory, accuracyDirectory, cmDirectory,
                    saveModelObject) {
  #require packages
  require(spectrolab)
  require(caret)
  require(dplyr)
  require(rlist)
  require(matrixStats)
  require(mlbench)
  require(MASS)
  require(parallel)
  
  #load spectra and convert to matrix and dataframe
  spec_raw = spectra
  spec_mat = as.matrix(spectra)
  spec_df1 = as.data.frame(spec_raw)
  
  #combine relavant metadata
  spec_df2 = as.data.frame(spec_mat)
  spec_df2 = cbind(spec_df2, spec_df1[className])
  colnames(spec_df2)[colnames(spec_df2) == className] <- className
  uniqueNames = unique(spec_df1[[className]])
  
  if (include_age == TRUE) {
    spec_df2$age = spec_df1$age
  }

  #create data partition: 70% of data for training, 30% for testing
  inTrain <- caret::createDataPartition(
    y = spec_df2[[className]],
    p = .7,
    list = FALSE
  )
  
  training <- spec_df2[inTrain,]
  testing <- spec_df2[-inTrain,]
  
  #tune model: 10-fold cross-validation repeated 3 times
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 10,
    sampling = resampling,
    repeats = 3)
  
  #Fit model. Note max iterations set to 100000 to allow model convergence
  plsFit <- train(
    as.formula(paste(className, "~ .")),
    data = training,
    maxit = 100000,
    method = "pls",
    trControl = ctrl,
    tuneLength = ncomp)
  
  if (resampling == 'down') {
    accFileName = paste(paste('accuracy', className, toString(iteration), sep = "_"),
                        "rds", sep = ".")
    saveRDS(plsFit$results$Accuracy, paste(accuracyDirectory, accFileName,
                                           sep = '/'))
  }
  
  #file name - change if you'd prefer a different file name
  fileName = paste(paste('pls', className, resampling, toString(iteration), sep = "_"),
                   "rds", sep = ".")
  
  if (saveModelObjec) {
    saveRDS(plsFit, paste(modelDirectory, fileName, sep = "/"))
  }
  
  if (resampling == 'up') {
    plsClasses = predict(plsFit, newdata = testing)
    cm = confusionMatrix(data = plsClasses, as.factor(testing[[className]]))
    cmFileName = paste(paste('cm', toString(iteration), sep="_"), "rds",
                       sep = ".")
    saveRDS(cm, paste(cmDirectory, cmFileName, sep = "/"))
  }
  
}

saveRDS(runPlsda, "revised pls/functions/runPlsda.rds")













