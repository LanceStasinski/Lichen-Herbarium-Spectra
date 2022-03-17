################################################################################
# Plot confusion matrix
################################################################################

plotConfusionMatrix = function(matrix, directory, fileName) {
  library(corrplot)
  library(naniar)
  
  #format matrix for plotting
  df = as.data.frame(matrix)
  df = df %>% replace_with_na_all(condition = ~.x == 0)
  mat = as.matrix(df)
  rownames(mat) = rownames(matrix)
  colnames(mat) = colnames(matrix)
  
  cols = colorRampPalette(c('white', '#fe9929'))
  
  jpeg(filename = paste(directory, fileName, sep = "/"),
       width = 12, height = 12, units = 'in', res = 1200)
  par(mfrow = c(1,1))
  corrplot::corrplot(mat,
                     cl.pos = 'n',
                     method = 'square',
                     tl.col = 'black',
                     cl.lim = c(0,1),
                     na.label = 'square',
                     na.label.col = 'white',
                     addCoef.col = '#542788',
                     number.digits = 2,
                     number.cex = .7,
                     col = cols(10))
  dev.off()
}

saveRDS(plotConfusionMatrix, "revised pls/functions/plotConfusionMatrix.rds")
