################################################################################
# Extract necessary data from model objects
################################################################################

getData = function(directory) {
  directory = 'revised pls/data/confusion-matrices'
  require(rlist)
  matrices = list.files(path = directory)
  
  cmList = list()
  accuracies = c()
  for (i in 1:length(matrices)) {
    cm = readRDS(paste(directory, matrices[i], sep = '/'))
    cmList = list.append(cmList, as.matrix(cm))
    accuracies = append(accuracies, cm$overall[1])
  }
  
  cmMean = t(Reduce("+", cmList)/length(matrices))
  cmMean = cmMean/rowSums(cmMean)
  
  # get standard deviations
  getSds = function(list){
    n = length(list); 	   
    rc = dim(list[[1]]); 	   
    ar1 = array(unlist(list), c(rc, n)); 	   
    round(apply(ar1, c(1, 2), sd), 2); 	         
  }
  
  cmSD = t(getSds(cmList))
  cmSD = cmSD/rowSums(cmSD)
  rownames(cmSD) = rownames(as.matrix(cmList[[1]]))
  colnames(cmSD) = colnames(as.matrix(cmList[[1]]))
  
  metrics = list(accuracies = accuracies, overallAccuracy = mean(accuracies),
                 cmMean = cmMean, cmSD = cmSD)
  
  return(metrics)
}

saveRDS(getData, "revised pls/functions/getData.rds")