################################################################################
# Vip plots
################################################################################

plotVip = function(modelDirectory, saveDirectory, baseFileName) {
  require(caret)
  require(rlist)
  models = list.files(path = modelDirectory)
  firstModel = readRDS(paste(modelDirectory, models[1], sep = "/"))
  importance = varImp(firstModel)$importance
  uniqueNames = colnames(importance)
  vipList = list()
  
  for (i in 1:length(uniqueNames)) {
    vipList = list.append(vipList, assign(uniqueNames[i],
                                          matrix(nrow = nrow(importance))))
  }
  
  for (i in 1:length(models)) {
    model = readRDS(paste(modelDirectory, models[i], sep = "/"))
    vip = varImp(model)$importance
    for (j in 1:length(uniqueNames)) {
      vipList[[j]] = cbind(vipList[[j]], vip[uniqueNames[j]])
    }
  }
  
  vip_to_bar = function(vip) {
    vip = vip[, -1]
    vip_mean = rowMeans(vip)
    sorted = sort(vip_mean)
    best = sorted[length(sorted)-5:length(sorted)]
    bm = as.data.frame(as.matrix(rev(best)))
    bm$color = '#d8b365' 
    worst = sorted[1:5]
    wm = as.data.frame(as.matrix(rev(worst)))
    wm$color = '#5ab4ac'
    m = rbind(bm, wm)
    barplot(rev(m[,1]), horiz = T, main = names(vip)[1],
            names.arg = rev(rownames(m)) , col = m$color)
  }
  
  fileNum = 1
  startIdx = 1
  endIdx = 12
  while (startIdx < length(vipList)) {
    fileName = paste(paste(baseFileName, fileNum, sep = "_"), "jpeg", sep = ".")
    jpeg(filename = paste(saveDirectory, fileName, sep = "/"),
         width = 6, height = 8, units = 'in', res = 1200)
    par(las = 2, mfrow = c(3,4))
    for (i in startIdx:endIdx) {
      if (is.null(vip[[i]])) next
      vip_to_bar(vipList[[i]])
    }
    dev.off()
    startIdx = endIdx + 1
    endIdx = endIdx + 12
    fileNum = fileNum + 1
  }
}

saveRDS(plotVip, "revised pls/functions/plotVip.rds")

plotVip(modelDirectory = 'revised pls/data/upsampling', saveDirectory = "revised pls/images/vip_plots",
        baseFileName = "testVip")