#Caret PLSDA
################################################################################
#Set up
################################################################################

library(corrplot)
library(matrixStats)
library(naniar)
library(spectrolab)

################################################################################
#run plsda
################################################################################
spectra = readRDS('spectra/lichen_spectra.rds')
#spectra = spectra[, seq(400, 2400, 100)]

classify = readRDS("functions/plsda.rds")

pls = classify(spectra = spectra, 
               className = "scientificName",
               ncomp = 4, 
               resampling = 'up',
               n_iteration = 25,
               include_age = T)

saveRDS(pls, 'models/morphology.rds')

pls = readRDS('models/species_age.rds')
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
pls = readRDS('models/class.rds')
a.fit = pls[[2]]
a.total = a.fit[,-1]
a.avg = as.matrix(rowMeans(a.total))
a.sd = as.matrix(rowSds(a.total))

a.lower = a.avg - a.sd
a.higher = a.avg + a.sd

lower.avgs = a.avg[a.avg < a.lower[which.max(a.avg)]]

#Graph to visually choose optimal number of components
jpeg(filename = '../../lichen figures/acc_vs_comp/class.jpeg',
     width = 10, height = 8, units = 'in', res = 600)
x = seq(1: length(a.avg))
par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
plot(x, a.avg, type = 'p', pch = 16, cex = .5, ylab = 'Accuracy', 
     xlab = 'Component', xlim = c(1,length(a.avg)), main = 'Class', 
     ylim = c(0,1))
arrows(x, a.lower, x, a.higher,length=0.05, angle=90, code=3)
abline(v = which.max(lower.avgs) + 1, col = 'blue')
abline(h = max(a.lower), col = "Red")
#legend('bottomright', legend = c('Mean', 'Within 1 SD of Max component','Optimal component'), 
       #pch = c(16, NA, NA), lty = c(NA, 1, 1), col = c('black', 'red', 'blue'))
dev.off()

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
write.csv(cm.sd, file = 'figures/confusion_matrices/standard deviations/class_age_sd.csv')

#format matrix for plotting
cm.total = as.data.frame(cm.total)
cm.total = cm.total %>% replace_with_na_all(condition = ~.x == 0)
cm.total = as.matrix(cm.total)
rownames(cm.total) = rownames(as.matrix(cm.list[[1]]))
colnames(cm.total) = colnames(as.matrix(cm.list[[1]]))

#save confusion matrix
write.csv(cm.total, "figures/confusion_matrices/cm_csv/class_age.csv")


#plot confusion matrix
cols = colorRampPalette(c('white', '#fe9929'))

jpeg(filename = '../../lichen figures/species-age_corrplot.jpeg',
     width = 12, height = 12, units = 'in', res = 1200)
par(mfrow = c(1,1))
corrplot::corrplot(as.matrix(cm.total),
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

mtext("Reference", side = 2, line = -8, cex = 2.5)
mtext("Prediction", side = 3, cex = 2.5, at = 2, line = 3)

################################################################################
#Variable importance
################################################################################
pls = readRDS('models/class_age.rds')

#plot vip as spectra
vip_to_spec = function(x){
  t.vip = t(x[,-1])
  t.vip = as.data.frame(t.vip)
  names(t.vip)[names(t.vip) == 'age'] = 2401
  colnames(t.vip) <- gsub("`", "", colnames(t.vip))
  s.vip = as_spectra(t.vip)
  plot(mean(s.vip), lwd = 1.5, lty = 1, ylim = c(0, 100),
       ylab = "Variable Importance", xlab = "Wavelength (nm)", cex.lab = 1,
       main = names(s.vip)[1])
  plot_quantile(s.vip, total_prob = 0.95, col = rgb(0, 0, 0, 0.25), 
                border = FALSE, add = TRUE)
}

jpeg(filename = '../../lichen figures/class_age_vip.jpeg',
     width = 12, height = 10, units = 'in', res = 1200)
vip.list = pls[[1]]
par(mfrow = c(2,3))
for (j in 1:length(vip.list)) {
  vip_to_spec(vip.list[[j]])
}
dev.off()

#plot top and bottom ten vip
vip.list = pls[[1]]
#vip = vip.list[[1]]
vip_to_bar = function(vip) {
  vip = vip[, -1]
  vip_mean = rowMeans(vip)
  sorted = sort(vip_mean)
  best = sorted[1998:2002]
  bm = as.data.frame(as.matrix(rev(best)))
  bm$color = '#d8b365' 
  worst = sorted[1:5]
  wm = as.data.frame(as.matrix(rev(worst)))
  wm$color = '#5ab4ac'
  m = rbind(bm, wm)
  barplot(rev(m[,1]), horiz = T, main = names(vip)[1],
          names.arg = rev(rownames(m)) , col = m$color)
}

jpeg(filename = '../../lichen figures/vip_rank_class.jpeg',
     width = 6, height = 8, units = 'in', res = 1200)
par(las=2, mfrow=c(2,3))
for (j in 1:19) {
  vip_to_bar(vip.list[[j]])
}
dev.off()
