################################################################################
#Caret PLSDA single function
##This function takes in spectra, the name of the column to be
#classified, the number of components to use, the type of resampling to be 
#done ('up' or 'down'), and the number of iterations to complete. The function 
#will return a list object that contains a list of matrices containing variable 
#importance values for each species, a matrix of model accuracy 
#(rows = components, columns = iteration), and vectors for overall accuracy and 
#kappa statistics. 
################################################################################

classify = function(spectra, className, ncomp, resampling, n_iteration, include_age) {
  #require packages
  require(spectrolab)
  require(caret)
  require(dplyr)
  require(rlist)
  require(matrixStats)
  require(mlbench)
  
  #load spectra and convert to matrix and dataframe
  spec_all = spectra
  spec_mat = as.matrix(spec_all)
  spec_all.df = as.data.frame(spec_all)
  

  #combine relevant meta data to matrix
  spec_df = as.data.frame(spec_mat)
  spec_df = cbind(spec_df, spec_all.df[className])
  colnames(spec_df)[colnames(spec_df) == className] <- className
  uniqueNames = unique(spec_all.df[[className]])
  
  if (include_age == TRUE) {
    spec_df$age = spec_all.df$age
  }

  ##################
  #Run PLSDA
  ##################
  
  #create vectors, lists, and matrices to store metrics and variable importance
  accuracy = c()
  kappa = c()
  a.fit = matrix(nrow = ncomp)
  cm.list = list()
  vip.list = list()
  results = list()
  model.list = list()

  #create variable importance matrix for each class
  for(j in 1:length(uniqueNames)){
    name = paste(uniqueNames[j], "vip", sep = ".")
    vip.list = list.append(vip.list, assign(name, matrix(nrow = ncol(spec_df)-1)))
  }
  

  
  #start of PLSDA code
  for(i in 1:n_iteration){
    
    #create data partition: 70% of data for training, 30% for testing
    inTrain <- caret::createDataPartition(
      y = spec_df[[className]],
      p = .7,
      list = FALSE
    )
    
    training <- spec_df[inTrain,]
    testing <- spec_df[-inTrain,]
    
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
  
    
    #variable importance
    if (length(uniqueNames) > 2) {
      vip = varImp(plsFit)
      for (k in 1:length(uniqueNames)) {
        class.vip = assign(paste0(uniqueNames[k], i), vip$importance[uniqueNames[k]])
        vip.list[[k]] = cbind(vip.list[[k]], get('class.vip'))
      }
    }
    
    
    #accuracy objects for determining n components
    a = assign(paste0('a', i), as.matrix(plsFit$results$Accuracy))
    a.fit <- cbind(a.fit, get('a'))
    
    #test model using the testing data partition (20% of data)
    plsClasses <- predict(plsFit, newdata = testing)
    
    #confusion/classification matrix objects to assess accuracy 
    cm = confusionMatrix(data = plsClasses, as.factor(testing[[className]]))
    cm.m = assign(paste0("cm", i), as.matrix(cm))
    cm.list <- list.append(cm.list, get('cm.m'))
    
    ac <- assign(paste0('acc',i), cm$overall[1])
    accuracy <- append(accuracy, get('ac'))
    
    kap = assign(paste0("kap",i), cm$overall[2])
    kappa <- append(kappa, get('kap'))
    

  }
  
  results = list.append(results, vip.list)
  results = list.append(results, a.fit)
  results = list.append(results, cm.list)
  results = list.append(results, accuracy)
  results = list.append(results, kappa)

  return(results)
  
}

saveRDS(classify, "functions/plsda.rds")

