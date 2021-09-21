#Linear Mixed Models
################################################################################
#setup
################################################################################

library(spectrolab)
library(lme4)
library(nlme)
library(ggplot2)
library(ggeffects)  
library(sjPlot)
library(rlist)
library(optimx)

################################################################################
#Evaluate optimizers
################################################################################
library(tidyverse)
library(lme4)
library(optimx)
library(parallel)
library(minqa)
library(dfoptim)

#Prepare data
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
data = meta(spectra)
spec.m = as.matrix(spectra) * 100
spectra_percent = as_spectra(spec.m)
meta(spectra_percent) = data
spec_df = as.data.frame(spectra_percent)

#model

lmm = lmer(spec_df[, '1932'] ~ age + (age|Class/Order/Family/scientificName), REML = T,
           data = spec_df)

#Run model with all optimizers
ncores <- detectCores() - 4
diff_optims <- allFit(lmm, maxfun = 1e5, parallel = 'multicore', ncpus = ncores)

#Output
is.OK <- sapply(diff_optims, is, "merMod")
diff_optims.OK <- diff_optims[is.OK]
lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)

################################################################################
#Compare models
################################################################################
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
data = meta(spectra)
spec.m = as.matrix(spectra) * 100
spectra_percent = as_spectra(spec.m)
meta(spectra_percent) = data
spec_df = as.data.frame(spectra_percent)

full_aic = c()
sound_aic = c()
species_aic = c()

full_bic = c()
sound_bic = c()
species_bic = c()

for(i in seq(400, 2400, 10)) {
    x = toString(i)
    full_hier = lmer(spec_df[, x] ~ age + (age|Class/Order/Family/scientificName),
                        data = spec_df, REML = T, 
                        lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
    full_aic = append(full_aic, AIC(full_hier))
    full_bic = append(full_bic, BIC(full_hier))
    
    sound_hier = lmer(spec_df[, x] ~ age + (age|Class/Order/Family),
                     data = spec_df, REML = T, 
                     lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
    sound_aic = append(sound_aic, AIC(sound_hier))
    sound_bic = append(sound_bic, BIC(sound_hier))
    
    species = lmer(spec_df[, x] ~ age + (age|scientificName),
                        data = spec_df, REML = T, 
                        lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
    species_aic = append(species_aic, AIC(species))
    species_bic = append(species_bic, BIC(species))
}

wv = seq(400, 2400, 10)
plot(wv, full_aic, col = 'blue', type = 'l', xlab = 'Wavelength (nm)', ylab = 'AIC')
lines(wv, sound_aic, col = 'red')
lines(wv, species_aic, col = 'black')
legend('bottomright',
       c('Full', 'Sound', 'Species'),
       col = c('blue', 'red', 'black'), lty = c(1,1))

wv = seq(400, 2400, 10)
plot(wv, full_bic, col = 'blue', type = 'l', xlab = 'Wavelength (nm)', ylab = 'BIC')
lines(wv, sound_bic, col = 'red')
lines(wv, species_bic, col = 'black')
legend('bottomright',
       c('Full', 'Sound', 'Species'),
       col = c('blue', 'red', 'black'), lty = c(1,1))

################################################################################
# Species as a random effect
################################################################################
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
data = meta(spectra)
spec.m = as.matrix(spectra) * 100

spectra_percent = as_spectra(spec.m)
meta(spectra_percent) = data

spec_df = as.data.frame(spectra_percent)

#spec_df$age = scale(spec_df$age, center = TRUE, scale = TRUE)

intercepts = as.data.frame(matrix(nrow = 29))
rownames(intercepts) = sort(unique(spec_df$scientificName))
slopes = as.data.frame(matrix(nrow = 29))
rownames(slopes) = sort(unique(spec_df$scientificName))

intercept_variance = c()
slope_variance = c()
resid_variance = c()

fixed_intercept_list = c()
fixed_slope_list = c()

intercept_2.5_list = c()
intercept_97.5_list = c()
slope_97.5_list = c()
slope_2.5_list = c()


for(i in seq(400, 2400, 1)) {
    x = toString(i)
    lmm = lmer(spec_df[, x] ~ age + (age|scientificName), data = spec_df, REML = T, 
               control = lmerControl(
                   optimizer ='optimx', optCtrl=list(method='nlminb')))
    
    coefs = coef(lmm)
    intercepts = cbind(intercepts, coefs[[1]][,1])
    slopes = cbind(slopes, coefs[[1]][,2])
    
    variances = as.data.frame(VarCorr(lmm))
    intercept_variance = append(intercept_variance, variances[1,4])
    slope_variance = append(slope_variance, variances[2,4])
    resid_variance = append(resid_variance, variances[4,4])
    
    fixed_intercept_list = append(fixed_intercept_list, as.numeric(fixef(lmm)[1]))
    fixed_slope_list = append(fixed_slope_list, as.numeric(fixef(lmm)[2]))
    
    ci = confint(lmm, method = 'Wald')
    intercept_2.5_list = append(intercept_2.5_list, ci[5])
    intercept_97.5_list = append(intercept_97.5_list, ci[11])
    slope_2.5_list = append(slope_2.5_list, ci[6])
    slope_97.5_list = append(slope_97.5_list, ci[12])
}

stats_list = list()

intercepts = intercepts[,-1]
colnames(intercepts) = seq(400,2400, 1)
slopes = slopes[,-1]
colnames(intercepts) = seq(400, 2400, 1)

stats_list = list.append(stats_list, intercepts) #1
stats_list = list.append(stats_list, slopes) #2
stats_list = list.append(stats_list, intercept_variance) #3
stats_list = list.append(stats_list, slope_variance) #4
stats_list = list.append(stats_list, resid_variance) #5
stats_list = list.append(stats_list, fixed_intercept_list) #6
stats_list = list.append(stats_list, fixed_slope_list) #7
stats_list = list.append(stats_list, intercept_2.5_list) #8
stats_list = list.append(stats_list, intercept_97.5_list) #9
stats_list = list.append(stats_list, slope_2.5_list) #10
stats_list = list.append(stats_list, slope_97.5_list) #11

saveRDS(stats_list, 'models/lmms/lmm_60yrs.rds')


################################################################################
#Plot
################################################################################

stats_list = readRDS('models/lmms/lmm_60yrs.rds')

par(mfrow = c(1,1))
#slopes
wv = seq(400, 2400, 1)
plot(wv, stats_list[[7]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (% reflectance/year)',
     ylim = c(min(stats_list[[2]]), max(stats_list[[2]])),
     main = 'Slopes')
polygon(c(wv, rev(wv)), c(stats_list[[10]], rev(stats_list[[11]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[2]])){
  lines(wv, stats_list[[2]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[7]])

#intercepts
wv = seq(400, 2400, 1)
plot(wv, stats_list[[6]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (% reflectance)',
     ylim = c(min(stats_list[[1]]), max(stats_list[[1]])),
     main = 'Intercepts')
polygon(c(wv, rev(wv)), c(stats_list[[8]], rev(stats_list[[9]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[1]])){
  lines(wv, stats_list[[1]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[6]])

#variances
wv = seq(400, 2400, 1)
plot(wv, stats_list[[3]],
     ylim = c(min(stats_list[[4]]), max(stats_list[[3]])),
     main = 'Variance',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[4]], col = 'red')
lines(wv, stats_list[[5]], col = 'gray')
legend('bottomright',
       c('Intercept', 'Slope', 'Residual'),
       col = c('blue', 'red', 'gray'), lty = c(1,1,1))

plot(wv, stats_list[[4]], main = 'Slope variance', ylab = 'variance',
     xlab = 'Wavelength (nm)', type = 'l')

################################################################################
#Hierarchical model
################################################################################
# setup data
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
data = meta(spectra)
spec.m = as.matrix(spectra) * 100
spectra_percent = as_spectra(spec.m)
meta(spectra_percent) = data
spec_df = as.data.frame(spectra_percent)

#slope and intercept values for each level
intercepts_species = as.data.frame(matrix(nrow = 29))
rownames(intercepts) = sort(unique(spec_df$scientificName))
slopes_species = as.data.frame(matrix(nrow = 29))
rownames(slopes) = sort(unique(spec_df$scientificName))

intercepts_family = as.data.frame(matrix(nrow = 19))
rownames(intercepts) = sort(unique(spec_df$Family))
slopes_family = as.data.frame(matrix(nrow = 19))
rownames(slopes) = sort(unique(spec_df$Family))

intercepts_order = as.data.frame(matrix(nrow = 16))
rownames(intercepts) = sort(unique(spec_df$Order))
slopes_order = as.data.frame(matrix(nrow = 16))
rownames(slopes) = sort(unique(spec_df$Order))

intercepts_class = as.data.frame(matrix(nrow = 6))
rownames(intercepts) = sort(unique(spec_df$Class))
slopes_class = as.data.frame(matrix(nrow = 6))
rownames(slopes) = sort(unique(spec_df$Class))

# Variance per level 
intercept_variance_species = c()
intercept_variance_family = c()
intercept_variance_order = c()
intercept_variance_class = c()

slope_variance_species = c()
slope_variance_family = c()
slope_variance_order = c()
slope_variance_class = c()

resid_variance = c()

#fixed slope and intercept values
fixed_intercept_list = c()
fixed_slope_list = c()

#95% confidence intervals for fixed effects
intercept_2.5_list = c()
intercept_97.5_list = c()
slope_97.5_list = c()
slope_2.5_list = c()


for(i in seq(400, 2400, 1)) {
  x = toString(i)
  lmm = lmer(spec_df[, '700'] ~ age + (age|Class/Order/Family/scientificName),
             data = spec_df, REML = T, 
             lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, 
                         optCtrl = list(maxfun = 1e5)))
  
  # coefficients for random effects of each level
  coefs = coef(lmm)
  intercepts_species = cbind(intercepts_species, coefs[[1]][,1])
  slopes_species = cbind(slopes_species, coefs[[1]][,2])
  
  intercepts_family = cbind(intercepts_family, coefs[[2]][,1])
  slopes_family = cbind(slopes_family, coefs[[2]][,2])
  
  intercepts_order = cbind(intercepts_order, coefs[[3]][,1])
  slopes_order = cbind(slopes_order, coefs[[3]][,2])
  
  intercepts_class = cbind(intercepts_class, coefs[[4]][,1])
  slopes_class = cbind(slopes_class, coefs[[4]][,2])
  
  #variance values per level
  variances = as.data.frame(VarCorr(lmm))
  
  intercept_variance_species = append(intercept_variance_species, variances[1,4])
  slope_variance_species = append(slope_variance_species, variances[2,4])
  
  intercept_variance_family = append(intercept_variance_family, variances[4,4])
  slope_variance_family = append(slope_variance_family, variances[5,4])
  
  intercept_variance_order = append(intercept_variance_order, variances[7,4])
  slope_variance_order = append(slope_variance_order, variances[8,4])
  
  intercept_variance_class = append(intercept_variance_class, variances[10,4])
  slope_variance_class = append(slope_variance_class, variances[11,4])
  
  resid_variance = append(resid_variance, variances[13,4])
  
  #fixed effects
  fixed_intercept_list = append(fixed_intercept_list, as.numeric(fixef(lmm)[1]))
  fixed_slope_list = append(fixed_slope_list, as.numeric(fixef(lmm)[2]))
  
  ci = confint(lmm, method = 'Wald')
  intercept_ci = ci['(Intercept)',]
  slope_ci = ci['age',]
  intercept_2.5_list = append(intercept_2.5_list, as.numeric(intercept_ci[1]))
  intercept_97.5_list = append(intercept_97.5_list, as.numeric(intercept_ci[2]))
  slope_2.5_list = append(slope_2.5_list, as.numeric(slope_ci[1]))
  slope_97.5_list = append(slope_97.5_list, as.numeric(slope_ci[2]))
}

stats_list = list()

clean = function(df) {
  df1 = df[, -1]
  colnames(df1) = seq(400, 2400, 1)
  return(df1)
}

intercepts_species = clean(intercepts_species)
intercepts_family = clean(intercepts_family)
intercepts_order = clean(intercepts_order)
intercepts_class = clean(intercepts_class)
slopes_species = clean(slopes_species)
slopes_family = clean(slopes_family)
slopes_order = clean(slopes_order)
slopes_class = clean(slopes_class)


stats_list = list.append(stats_list, intercepts_species) #1
stats_list = list.append(stats_list, intercepts_family) #2
stats_list = list.append(stats_list, intercepts_order) #3
stats_list = list.append(stats_list, intercepts_class) #4
stats_list = list.append(stats_list, slopes_species) #5
stats_list = list.append(stats_list, slopes_family) #6
stats_list = list.append(stats_list, slopes_order) #7
stats_list = list.append(stats_list, slopes_class) #8

stats_list = list.append(stats_list, intercept_variance_species) #9
stats_list = list.append(stats_list, intercept_variance_family) #10
stats_list = list.append(stats_list, intercept_variance_order) #11
stats_list = list.append(stats_list, intercept_variance_class) #12
stats_list = list.append(stats_list, slope_variance_species) #13
stats_list = list.append(stats_list, slope_variance_family) #14
stats_list = list.append(stats_list, slope_variance_order) #15
stats_list = list.append(stats_list, slope_variance_class) #16
stats_list = list.append(stats_list, resid_variance) #17

stats_list = list.append(stats_list, fixed_intercept_list) #18
stats_list = list.append(stats_list, fixed_slope_list) #19
stats_list = list.append(stats_list, intercept_2.5_list) #20
stats_list = list.append(stats_list, intercept_97.5_list) #21
stats_list = list.append(stats_list, slope_2.5_list) #22
stats_list = list.append(stats_list, slope_97.5_list) #23

saveRDS(stats_list, 'models/lmms/lmm_hier_60yrs.rds')

################################################################################
#Plot
################################################################################
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
data = meta(spectra)
spec.m = as.matrix(spectra) * 100
spectra_percent = as_spectra(spec.m)
meta(spectra_percent) = data
spec_df = as.data.frame(spectra_percent)
stats_list = readRDS('models/lmms/lmm_hier_60yrs.rds')

#species
par(mfrow = c(2,1))
##slopes
wv = seq(400, 2400, 1)
plot(wv, stats_list[[19]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (% reflectance/year)',
     ylim = c(min(stats_list[[5]]), max(stats_list[[5]])),
     main = 'Slopes - Species')
polygon(c(wv, rev(wv)), c(stats_list[[22]], rev(stats_list[[23]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[5]])){
  lines(wv, stats_list[[5]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[19]])

##intercepts
wv = seq(400, 2400, 1)
plot(wv, stats_list[[18]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (% reflectance)',
     ylim = c(min(stats_list[[1]]), max(stats_list[[21]])),
     main = 'Intercepts - Species')
polygon(c(wv, rev(wv)), c(stats_list[[20]], rev(stats_list[[21]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[1]])){
  lines(wv, stats_list[[1]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[18]])

################################################################################
#Family
par(mfrow = c(2,1))
##slopes
wv = seq(400, 2400, 1)
plot(wv, stats_list[[19]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (% reflectance/year)',
     ylim = c(min(stats_list[[6]]), max(stats_list[[23]])),
     main = 'Slopes - Family')
polygon(c(wv, rev(wv)), c(stats_list[[22]], rev(stats_list[[23]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[6]])){
  lines(wv, stats_list[[6]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[19]])

##intercepts
wv = seq(400, 2400, 1)
plot(wv, stats_list[[18]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (% reflectance)',
     ylim = c(min(stats_list[[2]]), max(stats_list[[2]])),
     main = 'Intercepts - Family')
polygon(c(wv, rev(wv)), c(stats_list[[20]], rev(stats_list[[21]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[2]])){
  lines(wv, stats_list[[2]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[18]])

################################################################################
#Order
par(mfrow = c(2,1))
##slopes
wv = seq(400, 2400, 1)
plot(wv, stats_list[[19]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (% reflectance/year)',
     ylim = c(min(stats_list[[22]]), max(stats_list[[23]])),
     main = 'Slopes - Order')
polygon(c(wv, rev(wv)), c(stats_list[[22]], rev(stats_list[[23]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[7]])){
  lines(wv, stats_list[[7]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[19]])

##intercepts
wv = seq(400, 2400, 1)
plot(wv, stats_list[[18]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (% reflectance)',
     ylim = c(min(stats_list[[3]]), max(stats_list[[21]])),
     main = 'Intercepts - Order')
polygon(c(wv, rev(wv)), c(stats_list[[20]], rev(stats_list[[21]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[3]])){
  lines(wv, stats_list[[3]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[18]])

################################################################################
#Class
par(mfrow = c(2,1))
##slopes
wv = seq(400, 2400, 1)
plot(wv, stats_list[[19]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (% reflectance/year)',
     ylim = c(min(stats_list[[8]]), max(stats_list[[8]])),
     main = 'Slopes - Class')
polygon(c(wv, rev(wv)), c(stats_list[[22]], rev(stats_list[[23]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[8]])){
  lines(wv, stats_list[[8]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[19]])

##intercepts
wv = seq(400, 2400, 1)
plot(wv, stats_list[[18]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (% reflectance)',
     ylim = c(min(stats_list[[4]]), max(stats_list[[4]])),
     main = 'Intercepts - Class')
polygon(c(wv, rev(wv)), c(stats_list[[20]], rev(stats_list[[21]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[4]])){
  lines(wv, stats_list[[4]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[18]])

################################################################################
#variances
par(mfrow = c(2,4))
##Species
wv = seq(400, 2400, 1)
plot(wv, stats_list[[9]],
     ylim = c(min(stats_list[[13]]), max(stats_list[[17]])),
     main = 'Variance - Species',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[13]], col = 'red')
lines(wv, stats_list[[17]], col = 'gray')

plot(wv, stats_list[[13]], main = 'Slope variance - Species', ylab = 'variance',
     xlab = 'Wavelength (nm)', type = 'l')

##Family
wv = seq(400, 2400, 1)
plot(wv, stats_list[[10]],
     ylim = c(min(stats_list[[14]]), max(stats_list[[10]])),
     main = 'Variance - Family',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[14]], col = 'red')
lines(wv, stats_list[[17]], col = 'gray')

plot(wv, stats_list[[14]], main = 'Slope variance - Family', ylab = 'variance',
     xlab = 'Wavelength (nm)', type = 'l')

##Order
wv = seq(400, 2400, 1)
plot(wv, stats_list[[11]],
     ylim = c(min(stats_list[[15]]), max(stats_list[[17]])),
     main = 'Variance - Order',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[15]], col = 'red')
lines(wv, stats_list[[17]], col = 'gray')

plot(wv, stats_list[[15]], main = 'Slope variance - Order', ylab = 'variance',
     xlab = 'Wavelength (nm)', type = 'l')

##Class
wv = seq(400, 2400, 1)
plot(wv, stats_list[[12]],
     ylim = c(min(stats_list[[16]]), max(stats_list[[12]])),
     main = 'Variance - Class',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[16]], col = 'red')
lines(wv, stats_list[[17]], col = 'gray')

plot(wv, stats_list[[16]], main = 'Slope variance - Class', ylab = 'variance',
     xlab = 'Wavelength (nm)', type = 'l')

##Total slope variance
par(mfrow = c(1,1))
wv = seq(400, 2400, 1)
plot(wv, stats_list[[13]],
     ylim = c(min(stats_list[[13]], stats_list[[14]], stats_list[[15]], stats_list[[16]]),
              max(stats_list[[13]], stats_list[[14]], stats_list[[15]], stats_list[[16]])),
     main = 'Variance - Slopes',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = '#a6cee3',
     type = 'l')
lines(wv, stats_list[[14]], col = '#1f78b4')
lines(wv, stats_list[[15]], col = '#b2df8a')
lines(wv, stats_list[[16]], col = '#33a02c')
legend('topright',
       c('Species', 'Family', 'Order', 'Class'),
       col = c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c'), lty = c(1,1,1,1))

##Total intercept variance
wv = seq(400, 2400, 1)
plot(wv, stats_list[[9]],
     ylim = c(min(stats_list[[9]], stats_list[[10]], stats_list[[11]], stats_list[[12]]),
              max(stats_list[[9]], stats_list[[10]], stats_list[[11]], stats_list[[12]])),
     main = 'Variance - Intercepts',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = '#a6cee3',
     type = 'l')
lines(wv, stats_list[[10]], col = '#1f78b4')
lines(wv, stats_list[[11]], col = '#b2df8a')
lines(wv, stats_list[[12]], col = '#33a02c')
legend('bottomright',
       c('Species', 'Family', 'Order', 'Class'),
       col = c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c'), lty = c(1,1,1,1))

##Total variance
wv = seq(400, 2400, 1)
plot(wv, stats_list[[9]],
     ylim = c(0, max(stats_list[[9]], stats_list[[10]], stats_list[[11]], stats_list[[12]])),
     main = 'Variance - Intercepts',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = '#fee090',
     type = 'l')
lines(wv, stats_list[[10]], col = '#fdae61')
lines(wv, stats_list[[11]], col = '#f46d43')
lines(wv, stats_list[[12]], col = '#d73027')
lines(wv, stats_list[[13]], col = '#e0f3f8')
lines(wv, stats_list[[14]], col = '#abd9e9')
lines(wv, stats_list[[15]], col = '#74add1')
lines(wv, stats_list[[16]], col = '#4575b4')
lines(wv, stats_list[[17]], col = 'gray')

legend('topright',
       c('Species - intercept', 'Family - intercept', 'Order - intercept', 'Class - intercept',
         'Species - slope', 'Family - slope', 'Order - slope', 'Class - slope', 'Residual'),
       col = c('#fee090', '#fdae61', '#f46d43', '#d73027',
               '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', 'gray'),
       lty = c(1,1,1,1,1,1,1,1,1))


################################################################################
#Hierarchical model 2
################################################################################
# setup data
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
data = meta(spectra)
spec.m = as.matrix(spectra) * 100
spectra_percent = as_spectra(spec.m)
meta(spectra_percent) = data
spec_df = as.data.frame(spectra_percent)

#slope and intercept values for each level


intercepts_family = as.data.frame(matrix(nrow = 19))
rownames(intercepts) = sort(unique(spec_df$Family))
slopes_family = as.data.frame(matrix(nrow = 19))
rownames(slopes) = sort(unique(spec_df$Family))

intercepts_order = as.data.frame(matrix(nrow = 16))
rownames(intercepts) = sort(unique(spec_df$Order))
slopes_order = as.data.frame(matrix(nrow = 16))
rownames(slopes) = sort(unique(spec_df$Order))

intercepts_class = as.data.frame(matrix(nrow = 6))
rownames(intercepts) = sort(unique(spec_df$Class))
slopes_class = as.data.frame(matrix(nrow = 6))
rownames(slopes) = sort(unique(spec_df$Class))

# Variance per level 
intercept_variance_family = c()
intercept_variance_order = c()
intercept_variance_class = c()

slope_variance_family = c()
slope_variance_order = c()
slope_variance_class = c()

resid_variance = c()

#fixed slope and intercept values
fixed_intercept_list = c()
fixed_slope_list = c()

#95% confidence intervals for fixed effects
intercept_2.5_list = c()
intercept_97.5_list = c()
slope_97.5_list = c()
slope_2.5_list = c()


for(i in seq(400, 2400, 10)) {
  x = toString(i)
  lmm = lmer(spec_df[, x] ~ age + (age|Class/Order/Family),
             data = spec_df, REML = T, 
             lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, 
                         optCtrl = list(maxfun = 1e5)))
  
  # coefficients for random effects of each level
  coefs = coef(lmm)
  
  intercepts_family = cbind(intercepts_family, coefs[[1]][,1])
  slopes_family = cbind(slopes_family, coefs[[1]][,2])
  
  intercepts_order = cbind(intercepts_order, coefs[[2]][,1])
  slopes_order = cbind(slopes_order, coefs[[2]][,2])
  
  intercepts_class = cbind(intercepts_class, coefs[[3]][,1])
  slopes_class = cbind(slopes_class, coefs[[3]][,2])
  
  #variance values per level
  variances = as.data.frame(VarCorr(lmm))
  
  
  intercept_variance_family = append(intercept_variance_family, variances[1,4])
  slope_variance_family = append(slope_variance_family, variances[2,4])
  
  intercept_variance_order = append(intercept_variance_order, variances[4,4])
  slope_variance_order = append(slope_variance_order, variances[5,4])
  
  intercept_variance_class = append(intercept_variance_class, variances[7,4])
  slope_variance_class = append(slope_variance_class, variances[8,4])
  
  resid_variance = append(resid_variance, variances[10,4])
  
  #fixed effects
  fixed_intercept_list = append(fixed_intercept_list, as.numeric(fixef(lmm)[1]))
  fixed_slope_list = append(fixed_slope_list, as.numeric(fixef(lmm)[2]))
  
  ci = confint(lmm, method = 'Wald')
  intercept_ci = ci['(Intercept)',]
  slope_ci = ci['age',]
  intercept_2.5_list = append(intercept_2.5_list, as.numeric(intercept_ci[1]))
  intercept_97.5_list = append(intercept_97.5_list, as.numeric(intercept_ci[2]))
  slope_2.5_list = append(slope_2.5_list, as.numeric(slope_ci[1]))
  slope_97.5_list = append(slope_97.5_list, as.numeric(slope_ci[2]))
}

stats_list = list()

clean = function(df) {
  df1 = df[, -1]
  colnames(df1) = seq(400, 2400, 10)
  return(df1)
}

intercepts_family = clean(intercepts_family)
intercepts_order = clean(intercepts_order)
intercepts_class = clean(intercepts_class)
slopes_family = clean(slopes_family)
slopes_order = clean(slopes_order)
slopes_class = clean(slopes_class)


stats_list = list.append(stats_list, intercepts_family) #1
stats_list = list.append(stats_list, intercepts_order) #2
stats_list = list.append(stats_list, intercepts_class) #3
stats_list = list.append(stats_list, slopes_family) #4
stats_list = list.append(stats_list, slopes_order) #5
stats_list = list.append(stats_list, slopes_class) #6

stats_list = list.append(stats_list, intercept_variance_family) #7
stats_list = list.append(stats_list, intercept_variance_order) #8
stats_list = list.append(stats_list, intercept_variance_class) #9
stats_list = list.append(stats_list, slope_variance_family) #10
stats_list = list.append(stats_list, slope_variance_order) #11
stats_list = list.append(stats_list, slope_variance_class) #12
stats_list = list.append(stats_list, resid_variance) #13

stats_list = list.append(stats_list, fixed_intercept_list) #14
stats_list = list.append(stats_list, fixed_slope_list) #15
stats_list = list.append(stats_list, intercept_2.5_list) #16
stats_list = list.append(stats_list, intercept_97.5_list) #17
stats_list = list.append(stats_list, slope_2.5_list) #18
stats_list = list.append(stats_list, slope_97.5_list) #19

saveRDS(stats_list, 'models/lmms/lmm_hier2_60yrs.rds')

################################################################################
#Plot
################################################################################
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
data = meta(spectra)
spec.m = as.matrix(spectra) * 100
spectra_percent = as_spectra(spec.m)
meta(spectra_percent) = data
spec_df = as.data.frame(spectra_percent)
stats_list = readRDS('models/lmms/lmm_hier_60yrs.rds')

#species
par(mfrow = c(2,1))
##slopes
wv = seq(400, 2400, 1)
plot(wv, stats_list[[19]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (% reflectance/year)',
     ylim = c(min(stats_list[[5]]), max(stats_list[[5]])),
     main = 'Slopes - Species')
polygon(c(wv, rev(wv)), c(stats_list[[22]], rev(stats_list[[23]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[5]])){
  lines(wv, stats_list[[5]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[19]])

##intercepts
wv = seq(400, 2400, 1)
plot(wv, stats_list[[18]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (% reflectance)',
     ylim = c(min(stats_list[[1]]), max(stats_list[[21]])),
     main = 'Intercepts - Species')
polygon(c(wv, rev(wv)), c(stats_list[[20]], rev(stats_list[[21]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[1]])){
  lines(wv, stats_list[[1]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[18]])

################################################################################
#Family
par(mfrow = c(2,1))
##slopes
wv = seq(400, 2400, 1)
plot(wv, stats_list[[19]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (% reflectance/year)',
     ylim = c(min(stats_list[[6]]), max(stats_list[[23]])),
     main = 'Slopes - Family')
polygon(c(wv, rev(wv)), c(stats_list[[22]], rev(stats_list[[23]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[6]])){
  lines(wv, stats_list[[6]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[19]])

##intercepts
wv = seq(400, 2400, 1)
plot(wv, stats_list[[18]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (% reflectance)',
     ylim = c(min(stats_list[[2]]), max(stats_list[[2]])),
     main = 'Intercepts - Family')
polygon(c(wv, rev(wv)), c(stats_list[[20]], rev(stats_list[[21]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[2]])){
  lines(wv, stats_list[[2]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[18]])

################################################################################
#Order
par(mfrow = c(2,1))
##slopes
wv = seq(400, 2400, 1)
plot(wv, stats_list[[19]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (% reflectance/year)',
     ylim = c(min(stats_list[[22]]), max(stats_list[[23]])),
     main = 'Slopes - Order')
polygon(c(wv, rev(wv)), c(stats_list[[22]], rev(stats_list[[23]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[7]])){
  lines(wv, stats_list[[7]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[19]])

##intercepts
wv = seq(400, 2400, 1)
plot(wv, stats_list[[18]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (% reflectance)',
     ylim = c(min(stats_list[[3]]), max(stats_list[[21]])),
     main = 'Intercepts - Order')
polygon(c(wv, rev(wv)), c(stats_list[[20]], rev(stats_list[[21]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[3]])){
  lines(wv, stats_list[[3]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[18]])

################################################################################
#Class
par(mfrow = c(2,1))
##slopes
wv = seq(400, 2400, 1)
plot(wv, stats_list[[19]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (% reflectance/year)',
     ylim = c(min(stats_list[[8]]), max(stats_list[[8]])),
     main = 'Slopes - Class')
polygon(c(wv, rev(wv)), c(stats_list[[22]], rev(stats_list[[23]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[8]])){
  lines(wv, stats_list[[8]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[19]])

##intercepts
wv = seq(400, 2400, 1)
plot(wv, stats_list[[18]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (% reflectance)',
     ylim = c(min(stats_list[[4]]), max(stats_list[[4]])),
     main = 'Intercepts - Class')
polygon(c(wv, rev(wv)), c(stats_list[[20]], rev(stats_list[[21]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[4]])){
  lines(wv, stats_list[[4]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[18]])

################################################################################
#variances
par(mfrow = c(2,4))
##Species
wv = seq(400, 2400, 1)
plot(wv, stats_list[[9]],
     ylim = c(min(stats_list[[13]]), max(stats_list[[17]])),
     main = 'Variance - Species',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[13]], col = 'red')
lines(wv, stats_list[[17]], col = 'gray')

plot(wv, stats_list[[13]], main = 'Slope variance - Species', ylab = 'variance',
     xlab = 'Wavelength (nm)', type = 'l')

##Family
wv = seq(400, 2400, 1)
plot(wv, stats_list[[10]],
     ylim = c(min(stats_list[[14]]), max(stats_list[[10]])),
     main = 'Variance - Family',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[14]], col = 'red')
lines(wv, stats_list[[17]], col = 'gray')

plot(wv, stats_list[[14]], main = 'Slope variance - Family', ylab = 'variance',
     xlab = 'Wavelength (nm)', type = 'l')

##Order
wv = seq(400, 2400, 1)
plot(wv, stats_list[[11]],
     ylim = c(min(stats_list[[15]]), max(stats_list[[17]])),
     main = 'Variance - Order',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[15]], col = 'red')
lines(wv, stats_list[[17]], col = 'gray')

plot(wv, stats_list[[15]], main = 'Slope variance - Order', ylab = 'variance',
     xlab = 'Wavelength (nm)', type = 'l')

##Class
wv = seq(400, 2400, 1)
plot(wv, stats_list[[12]],
     ylim = c(min(stats_list[[16]]), max(stats_list[[12]])),
     main = 'Variance - Class',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[16]], col = 'red')
lines(wv, stats_list[[17]], col = 'gray')

plot(wv, stats_list[[16]], main = 'Slope variance - Class', ylab = 'variance',
     xlab = 'Wavelength (nm)', type = 'l')

##Total slope variance
par(mfrow = c(1,1))
wv = seq(400, 2400, 1)
plot(wv, stats_list[[13]],
     ylim = c(min(stats_list[[13]], stats_list[[14]], stats_list[[15]], stats_list[[16]]),
              max(stats_list[[13]], stats_list[[14]], stats_list[[15]], stats_list[[16]])),
     main = 'Variance - Slopes',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = '#a6cee3',
     type = 'l')
lines(wv, stats_list[[14]], col = '#1f78b4')
lines(wv, stats_list[[15]], col = '#b2df8a')
lines(wv, stats_list[[16]], col = '#33a02c')
legend('topright',
       c('Species', 'Family', 'Order', 'Class'),
       col = c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c'), lty = c(1,1,1,1))

##Total intercept variance
wv = seq(400, 2400, 1)
plot(wv, stats_list[[9]],
     ylim = c(min(stats_list[[9]], stats_list[[10]], stats_list[[11]], stats_list[[12]]),
              max(stats_list[[9]], stats_list[[10]], stats_list[[11]], stats_list[[12]])),
     main = 'Variance - Intercepts',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = '#a6cee3',
     type = 'l')
lines(wv, stats_list[[10]], col = '#1f78b4')
lines(wv, stats_list[[11]], col = '#b2df8a')
lines(wv, stats_list[[12]], col = '#33a02c')
legend('bottomright',
       c('Species', 'Family', 'Order', 'Class'),
       col = c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c'), lty = c(1,1,1,1))

##Total variance
wv = seq(400, 2400, 1)
plot(wv, stats_list[[9]],
     ylim = c(0, max(stats_list[[9]], stats_list[[10]], stats_list[[11]], stats_list[[12]])),
     main = 'Variance - Intercepts',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = '#fee090',
     type = 'l')
lines(wv, stats_list[[10]], col = '#fdae61')
lines(wv, stats_list[[11]], col = '#f46d43')
lines(wv, stats_list[[12]], col = '#d73027')
lines(wv, stats_list[[13]], col = '#e0f3f8')
lines(wv, stats_list[[14]], col = '#abd9e9')
lines(wv, stats_list[[15]], col = '#74add1')
lines(wv, stats_list[[16]], col = '#4575b4')
lines(wv, stats_list[[17]], col = 'gray')

legend('topright',
       c('Species - intercept', 'Family - intercept', 'Order - intercept', 'Class - intercept',
         'Species - slope', 'Family - slope', 'Order - slope', 'Class - slope', 'Residual'),
       col = c('#fee090', '#fdae61', '#f46d43', '#d73027',
               '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', 'gray'),
       lty = c(1,1,1,1,1,1,1,1,1))

