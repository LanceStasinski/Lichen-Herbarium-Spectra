#Caret PLSDA
################################################################################
#Set up
################################################################################

library(corrplot)
library(matrixStats)
library(naniar)

################################################################################
#run plsda
################################################################################

classify = readRDS("functions/plsda.rds")

pls = classify(file = "spectra/lichen_spectra.rds", 
               className = "scientificName",
               ncomp = 3, 
               resampling = 'down',
               n_iteration = 3)

################################################################################
#Assess accuracy and kappa
################################################################################
accuracy = pls[[4]]
mean(accuracy)
sd(accuracy)

kappa = pls[[5]]
mean(kappa)
sd(kappa)

################################################################################
#Accuracy values for choosing the optimal number of components to use
################################################################################
a.fit = pls[[2]]
a.total = a.fit[,-1]
a.avg = as.matrix(rowMeans(a.total))
a.sd = as.matrix(rowSds(a.total))

a.lower = a.avg - a.sd
a.higher = a.avg + a.sd

#Graph to visually choose optimal number of components
x = 1:3
par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
plot(x, a.avg, type = 'p', pch = 16, cex = .75, ylab = 'Accuracy', 
     xlab = 'Component', xlim = c(1,60), main = 'Accuracy for Species_ID', 
     ylim = c(0,1))
arrows(x, a.lower, x, a.higher,length=0.05, angle=90, code=3)
abline(v = which.max(a.avg), col = 'blue')
abline(h = max(a.avg), col = "Red")
legend('bottomright', legend = c('Mean', 'Maximum accuracy','Best component'), 
       pch = c(16, NA, NA), lty = c(NA, 1, 1), col = c('black', 'red', 'blue'))

################################################################################
#Confusion/Classification Matrices
################################################################################
#take average of confusion matrices, reorient matrix, change averages to 
#proportions
cm.list = pls[[3]]
cm.avg = Reduce('+', cm.list)/100
cm.avg = t(cm.avg)
cm.total = cm.avg/rowSums(cm.avg)

#standard deviations
f1 <- function(lst){
  n <- length(lst); 	   
  rc <- dim(lst[[1]]); 	   
  ar1 <- array(unlist(lst), c(rc, n)); 	   
  round(apply(ar1, c(1, 2), sd), 2); 	         
}
cm.sd = f1(cm.list)
cm.sd = t(cm.sd)
cm.sd = cm.sd/rowSums(cm.avg)
rownames(cm.sd) = rownames(as.matrix(cm.list[[1]]))
colnames(cm.sd) = colnames(as.matrix(cm.list[[1]]))
write.csv(cm.sd, file = 'figures/confusion_matrices/standard deviations/Species_sd.csv')

#format matrix for plotting
cm.total = as.data.frame(cm.total)
cm.total = cm.total %>% replace_with_na_all(condition = ~.x == 0)
cm.total = as.matrix(cm.total)
rownames(cm.total) = rownames(as.matrix(cm.list[[1]]))
colnames(cm.total) = colnames(as.matrix(cm.list[[1]]))

#save confusion matrix
write.csv(cm.total, "figures/confusion_matrices/cm_csv/Species.csv")


#plot confusion matrix
cols = colorRampPalette(c('#f5f5f5', '#fe9929'))

par(mfrow = c(1,1))
corrplot::corrplot(as.matrix(cm.total),
                   is.corr = F, 
                   method = 'square', 
                   col = cols(10),
                   addCoef.col = '#542788',
                   tl.srt = 0, 
                   tl.offset = 1, 
                   number.digits = NULL, 
                   tl.cex = 1.2, 
                   cl.cex = 1, 
                   number.cex = 1.5,
                   tl.col = 'black', 
                   tl.pos = 'n',
                   cl.pos = 'n',
                   cl.lim = c(0,1),
                   na.label = 'square', 
                   na.label.col = 'white',
                   addgrid.col = 'grey')
mtext("Reference", side = 2, line = -8, cex = 2.5)
mtext("Prediction", side = 3, cex = 2.5, at = 2, line = 3)

################################################################################
#Variable importance
################################################################################
vip_to_spec = function(x){
  t.vip = t(x)
  colnames(t.vip) <- seq(400,2400, by = 1)
  s.vip = as_spectra(t.vip)
}
do.vip = do.vip[,-1]
do.vip.spec = vip_to_spec(do.vip)

da.vip = da.vip[,-1]
da.vip.spec = vip_to_spec(da.vip)

#dx.vip = dx.vip[,-1]
#dx.vip.spec = vip_to_spec(dx.vip)

#plot
par(mfrow = c(2,1))

plot(mean(da.vip.spec), lwd = 1.5, lty = 1, col = '#00B0F6', ylim = c(0, 100),
     ylab = "Variable Importance", xlab = NA, cex.lab = 1.5)
plot_quantile(da.vip.spec, total_prob = 0.95, col = rgb(0, 0.69, 0.965, 0.25), 
              border = FALSE, add = TRUE)


plot(mean(do.vip.spec), lwd = 2, lty = 1, col = '#F8766D', 
     cex.lab = 1.5, ylim = c(0, 100), ylab = "Variable Importance", 
     xlab = NA)
plot_quantile(do.vip.spec, total_prob = 0.95, col = rgb(0.972549,0.4627451,0.427451, 0.25), 
              border = FALSE, add = TRUE)


plot(mean(dx.vip.spec), lwd = 1.5, lty = 1, col = rgb(0,0,0,1), ylim = c(0, 100),
     ylab = "Variable Importance", xlab = 'Wavelength (nm)', cex.lab = 1.5)
plot_quantile(dx.vip.spec, total_prob = 0.95, col = rgb(0, 0, 0, 0.25), 
              border = FALSE, add = TRUE)
