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
#Report desired statistics from all wavelengths
################################################################################
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
data = meta(spectra)
spec.m = as.matrix(spectra) * 100

spectra_percent = as_spectra(spec.m)
meta(spectra_percent) = data

spec_df = as.data.frame(spectra_percent)
#spec_df$age = scale(spec_df$age, center = TRUE, scale = TRUE)

varSlope_aic = c()
fixedSlope_aic = c()

for(i in seq(400, 2400, 10)) {
    x = toString(i)
    lmm_varSlope = lmer(spec_df[, x] ~ age + (age|Class/Order/Family/scientificName),
                        data = spec_df, REML = T, 
                        lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
    lmm_fixedSlope = lmer(spec_df[, x] ~ age + (1|Class/Order/Family/scientificName),
                          data = spec_df, REML = T, 
                          lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
    
    varSlope_aic = append(varSlope_aic, AIC(lmm_varSlope))
    fixedSlope_aic = append(fixedSlope_aic, AIC(lmm_fixedSlope))
}

wv = seq(400, 2400, 10)
plot(wv, varSlope_aic, col = 'blue', type = 'l', xlab = 'Wavelength (nm)', ylab = 'AIC')
lines(wv, fixedSlope_aic, col = 'red')
legend('bottomright',
       c('Variable slope & intercept', 'Variable intercept'),
       col = c('blue', 'red'), lty = c(1,1))

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
  lmm = lmer(spec_df[, x] ~ age + (age|Class/Order/Family/scientificName),
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
  intercept_2.5_list = append(intercept_2.5_list, ci[5])
  intercept_97.5_list = append(intercept_97.5_list, ci[11])
  slope_2.5_list = append(slope_2.5_list, ci[6])
  slope_97.5_list = append(slope_97.5_list, ci[12])
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
stats_list = list.append(stats_list, resid_variance) #5

stats_list = list.append(stats_list, fixed_intercept_list) #6
stats_list = list.append(stats_list, fixed_slope_list) #7
stats_list = list.append(stats_list, intercept_2.5_list) #8
stats_list = list.append(stats_list, intercept_97.5_list) #9
stats_list = list.append(stats_list, slope_2.5_list) #10
stats_list = list.append(stats_list, slope_97.5_list) #11

saveRDS(stats_list, 'models/lmms/lmm_hier_60yrs.rds')





