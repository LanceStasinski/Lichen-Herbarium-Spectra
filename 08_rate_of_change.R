#Rate of change
################################################################################
#setup
################################################################################

library(spectrolab)
library(rlist)
library(splines)
################################################################################
#comparison of functions - testing
################################################################################

spectra = readRDS('spectra/lichen_spectra.rds')

spec = spectra[meta(spectra)$scientificName == 'Flavoparmelia_caperata',]
spec700 = spec[, 700]
spec_df = as.data.frame(spec)
spec_m = as.data.frame(as.matrix(spec700))
colnames(spec_m) =  'ref'
spec_m$age = spec_df$age
spec_m = spec_m[spec_m$age <= 60,]

linear = lm(ref~age, data = spec_m)
expon = lm(log(ref)~age, data = spec_m)
sp = lm(ref ~ ns(age, df = 6), data = spec_m)


timevalues <- seq(0, 123, 1)
expon2 <- exp(predict(expon, list(age=timevalues)))
linear2 = predict(linear, list(age=timevalues))

plot(spec_m$age, spec_m$ref, pch=16, xlab = "Age (years)", ylab = "Reflectance")
lines(timevalues, expon2,lwd=2, col = "red")
lines(timevalues, linear2, lwd=2, col = 'blue')
lines(timevalues, predict(sp, data.frame(age = timevalues)), col = 'green')

AIC(linear)
AIC(expon)
AIC(sp)
################################################################################
#comparison of functions - all samples
################################################################################
spectra = readRDS('spectra/lichen_spectra.rds')
spec_df = as.data.frame(spectra)


lin_aic = c()
ex_aic = c()
sp_aic = c()

for(i in seq(400, 2400, 1)) {
  x = toString(i)
  linear = lm(spec_df[, x] ~ spec_df$age)
  expon = lm(log(spec_df[, x]) ~ spec_df$age)
  sp = lm(spec_df[, x] ~ ns(spec_df$age, df = 6))
  
  lin_aic = append(lin_aic, AIC(linear))
  ex_aic = append(ex_aic, AIC(expon))
  sp_aic = append(sp_aic, AIC(sp))
}


max = max(c(max(lin_aic), max(ex_aic), max(sp_aic)))
min = min(c(min(lin_aic), min(ex_aic), min(sp_aic)))

wv = seq(400, 2400, 1)

#plot all: linear, exponential, spline
plot(wv, lin_aic, type='l', lty = 1, col = 'blue', ylim = c(min, max), ylab = 'AIC', xlab = 'Wavelength (nm)')
lines(wv, ex_aic, lty = 1, col = 'red' )
lines(wv, sp_aic, lty = 1, col = 'green')
legend('topright',
       legend = c('Linear', 'Exponential', 'Spline'),
       lty = 1, col = c('blue', 'red', 'green'))

#plot linear and spline AIC
plot(wv, lin_aic, type='l', lty = 1, col = 'blue', ylab = 'AIC', xlab = 'Wavelength (nm)')
lines(wv, sp_aic, lty = 1, col = 'green')
legend('topright',
       legend = c('Linear', 'Spline'),
       lty = 1, col = c('blue', 'green'))


################################################################################
#comparison of functions - by taxonomic unit
################################################################################
functionTypeComparison = function(spectra, splineDF, taxa){
  
  spec_df = as.data.frame(spectra)
  uniqueNames = unique(spec_df[, taxa])
  
  completeList = list()
  nameList = c()
  fullAicList = list()
  
  for (i in 1:length(uniqueNames)) {
    species_spec = spec_df[spec_df[, taxa] == uniqueNames[i],]
    
    if (length(unique(species_spec$age)) < splineDF + 2) {
      next
    }
    
    nameList = append(nameList, uniqueNames[i])
    
    lin_aic = c()
    ex_aic = c()
    sp_aic = c()
    
    for(j in seq(400, 2400, 1)) {
      wl = toString(j)
      linear = lm(species_spec[, wl] ~ species_spec$age)
      expon = lm(log(species_spec[, wl]) ~ species_spec$age)
      sp = lm(species_spec[, wl] ~ ns(species_spec$age, df = splineDF))
      
      lin_aic = append(lin_aic, AIC(linear))
      ex_aic = append(ex_aic, AIC(expon))
      sp_aic = append(sp_aic, AIC(sp))
    }
    
    taxonAIC = Reduce(cbind, list(lin_aic, ex_aic, sp_aic))
    colnames(taxonAIC) = c('linear', 'expon', 'spline')
    
    fullAicList = list.append(fullAicList, taxonAIC)
    
  }
  
  completeList = list.append(completeList, fullAicList)
  completeList = list.append(completeList, nameList)
  
  return(completeList)
}

spectra = readRDS('spectra/lichen_spectra.rds')
taxaAIC = functionTypeComparison(spectra = spectra, splineDF = 6, taxa = 'scientificName')

spec_df = as.data.frame(spectra)
taxa = 'scientificName'

par(mfrow = c(3,5))
for(x in 1:length(taxaAIC[[2]])) {
  taxon = as.data.frame(taxaAIC[[1]][x])
  max = max(c(max(taxon$linear), max(taxon$expon), max(taxon$spline)))
  min = min(c(min(taxon$linear), min(taxon$expon), min(taxon$spline)))
  wv = seq(400, 2400, 1)
  
  #plot all: linear, exponential, spline
  #plot(wv, taxon$linear, type='l', lty = 1, col = 'blue', ylim = c(min, max),
       #ylab = 'AIC', xlab = 'Wavelength (nm)', main = taxaAIC[[2]][x])
  #lines(wv, taxon$expon, lty = 1, col = 'red' )
  #lines(wv, taxon$spline, lty = 1, col = 'green')
  
  #plot linear and spline AIC
  plot(wv, taxon$linear, type='l', lty = 1, col = 'blue', ylab = 'AIC', xlab = 'Wavelength (nm)', main = taxaAIC[[2]][x])
  lines(wv, taxon$spline, lty = 1, col = 'green')
  
}


################################################################################
#plot - may delete soon
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
sp = lm(reflectance ~ ns(age, df = 6), data = spec)
lines(x, predict(sp, data.frame(age = x)), col = 'blue')

v = plot(spec$age, spec$reflectance)
abline(linear, col = 'red')
x = seq(0, 120, 1)
lines(x, predict(logit, data.frame(age = x), type="response"), col = 'blue')
lines(x, predict(sp, data.frame(age = x)), col = 'green')

################################################################################
#slope by age classes
################################################################################ 

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

spectra = readRDS('spectra/lichen_spectra.rds')
spectra = normalize(spectra)
spectra = spectra[meta(spectra)$age < 60,]

uniqueNames = unique(meta(spectra)$Class)
speciesList = list()
nameList = c()

for(i in 1:length(uniqueNames)) {
  classSpec = spectra[meta(spectra)$Class == uniqueNames[i],]
  classSpec_df = as.data.frame(classSpec)

  #skip if time series has less than 12 specimens (12 specimens * 4 scans = 48)
  if(nrow(classSpec_df) < 48) {
    next
  }
  nameList = append(nameList, unique(classSpec_df$Class))
  
  df.list = list()
  #break age into 10 year intervals
  for (j in seq(20, 60, 20)) {
    ageSpec = classSpec[meta(classSpec)$age %in% seq(j-19, j, 1)] 
    ageSpec_df = as.data.frame(ageSpec)
    #calculate slope for each wavelength
    slope.list = c()
    p.list = c()
    for (y in 400:2400) {
      spec = as.data.frame(ageSpec[,y])
      spec = cbind(spec,ageSpec_df$age)
      colnames(spec) = c('reflectance', 'age')
      m = lm(reflectance ~ age, data = spec)
      slope.list = append(slope.list, m$coefficients[2])
      p.list = append(p.list, lmp(m))
    }
  
  
  slope = as.data.frame(slope.list)
  colnames(slope) = c('reflectance')
  slope$wavelength = seq(400, 2400, by = 1)
  slope = cbind(slope, p.list)
  slope$color = 'black'
  slope$color[slope$p.list > 0.05] = 'gray'
  #append df to list of dfs for this species
  df.list = list.append(df.list, slope)
  }
  speciesList = list.append(speciesList, df.list)
}

#plot
par(mfrow = c(1, 3))
for(x in 1:length(speciesList)) {
  species = speciesList[[x]]
  
  maxes = c()
  mins = c()
  for (z in 1:length(species)) {
    max = max(species[[z]]$reflectance)
    maxes = append(maxes, max)
    min = min(species[[z]]$reflectance)
    mins = append(mins, min)
  }
  
  plot(species[[1]]$wavelength,
       species[[1]]$reflectance,
       type = 'l',
       lty = 1,
       lwd = 2,
       col = species[[1]]$color,
       main = nameList[x],
       xlab = 'Wavelength (nm)',
       ylab = 'Change in reflectance per year',
       ylim = c(min(mins), max(maxes)))
  lines(species[[2]]$wavelength,
        species[[2]]$reflectance,
        col = species[[2]]$color,
        lty = 2,
        lwd = 2)
  lines(species[[3]]$wavelength,
        species[[3]]$reflectance,
        col = species[[3]]$color,
        lty = 3,
        lwd = 2)
  abline(h=0, lty = 6, col = 'blue')
  
}


################################################################################
#Plot all samples by age class
################################################################################

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

spectra = readRDS('spectra/lichen_spectra.rds')
spectra = normalize(spectra)
spectra = spectra[meta(spectra)$age < 60,]


classSpec = spectra
classSpec_df = as.data.frame(classSpec)


df.list = list()
#break age into 10 year intervals
for (j in seq(20, 60, 20)) {
  ageSpec = classSpec[meta(classSpec)$age %in% seq(j-19, j, 1)] 
  ageSpec_df = as.data.frame(ageSpec)
  #calculate slope for each wavelength
  slope.list = c()
  p.list = c()
  for (y in 400:2400) {
    spec = as.data.frame(ageSpec[,y])
    spec = cbind(spec,ageSpec_df$age)
    colnames(spec) = c('reflectance', 'age')
    m = lm(reflectance ~ age, data = spec)
    slope.list = append(slope.list, m$coefficients[2])
    p.list = append(p.list, lmp(m))
  }
  
  
  slope = as.data.frame(slope.list)
  colnames(slope) = c('reflectance')
  slope$wavelength = seq(400, 2400, by = 1)
  slope = cbind(slope, p.list)
  slope$color = 'black'
  slope$color[slope$p.list > 0.05] = 'gray'
  #append df to list of dfs for this species
  df.list = list.append(df.list, slope)
}


#plot
par(mfrow = c(1, 1))

species = df.list

maxes = c()
mins = c()
for (z in 1:length(species)) {
  max = max(species[[z]]$reflectance)
  maxes = append(maxes, max)
  min = min(species[[z]]$reflectance)
  mins = append(mins, min)
}

plot(species[[1]]$wavelength,
     species[[1]]$reflectance,
     type = 'l',
     lty = 1,
     lwd = 2,
     col = species[[1]]$color,
     main = 'All samples (normalized)',
     xlab = 'Wavelength (nm)',
     ylim = c(min(mins), max(maxes)),
     ylab = 'Change in reflectance per year')
lines(species[[2]]$wavelength,
      species[[2]]$reflectance,
      col = species[[2]]$color,
      lty = 2,
      lwd = 2)
lines(species[[3]]$wavelength,
      species[[3]]$reflectance,
      col = species[[3]]$color,
      lty = 3,
      lwd = 2)
abline(h=0, lty = 6, col = 'blue')
legend('topright',
       legend = c('1-20 years', '21-40 years', '41-60 years'),
       lty = 1:3)


