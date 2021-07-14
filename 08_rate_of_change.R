#Rate of change
################################################################################
#setup
################################################################################

library(spectrolab)
library(rlist)

################################################################################
#rate of per wavelength per class
################################################################################

spectra = readRDS('spectra/lichen_spectra.rds')

uniqueNames = unique(meta(spectra)$scientificName)
classList = list()

for(i in 1:length(uniqueNames)) {
  classSpec = spectra[meta(spectra)$scientificName == 'Flavoparmelia_caperata',]#uniqueNames[1],]
  classSpec_df = as.data.frame(classSpec)
  slope.list = c()
  for (j in 400:2400) {
    spec = as.data.frame(classSpec[,400])
    spec = cbind(spec, classSpec_df$age)
    colnames(spec) = c('reflectance', 'age')
    m = lm(reflectance ~ age, data = spec)
    slope.list = append(slope.list, m$coefficients[2])
  }
  
  slope = as.data.frame(slope.list)
  colnames(slope) = c('reflectance')
  slope$wavelength = seq(400, 2400, by = 1)
  classList = list.append(classList, slope)
}

saveRDS(classList, 'models/slopes/species.rds')

################################################################################
#plot
################################################################################
par(mfrow = c(3, 5))
for (i in 1:length(uniqueNames)) {
  if (is.na(classList[[i]]$reflectance[[1]] == TRUE)) {
    next
  }
  plot(classList[[i]]$wavelength, classList[[i]]$reflectance,
       ylab = 'Change in reflectance per year',
       xlab = 'Wavelength (nm)',
       main = uniqueNames[i],
       ylim = c(-.015, .05))
       
}


#models
spectra = readRDS('spectra/lichen_spectra.rds')
uniqueNames = unique(meta(spectra)$scientificName)
classList = list()
classSpec = spectra[meta(spectra)$scientificName == 'Flavoparmelia_caperata',]#uniqueNames[1],]
classSpec_df = as.data.frame(classSpec)
spec = as.data.frame(classSpec[,400])
spec = cbind(spec, classSpec_df$age)
colnames(spec) = c('reflectance', 'age')

linear = lm(reflectance ~ age, data = spec)
logit = glm(reflectance ~ age, data = spec, family = binomial)
sp = lm(reflectance ~ ns(age, df = 3), data = spec)
lines(x, predict(sp, data.frame(age = x)), col = 'blue')
