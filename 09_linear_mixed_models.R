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
library(forecast)
library(partR2)


#example
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
spectra = aggregate(spectra, by = meta(spectra)$X, mean, try_keep_txt(mean))
data = meta(spectra)
spec.m = as.matrix(spectra) * 100
spectra_percent = as_spectra(spec.m)
meta(spectra_percent) = data
spec_df = as.data.frame(spectra_percent)

lmm1 = lmer(spec_df[, '700'], ~ age + (age|Class/Order/Family/scientificName),
           data = spec_df, REML = T, 
           lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, 
                       optCtrl = list(maxfun = 1e5)))
summary(lmm1)
plot(lmm1)

lmm2 = lmer(spec_df[, '700'] ~ age + (age|Class/Order/Family),
           data = spec_df, REML = T, 
           lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, 
                       optCtrl = list(maxfun = 1e5)))
summary(lmm2)
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

lmm = lmer(spec_df[, '700'] ~ age + (age|Class/Order/Family/scientificName), REML = T,
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
#spectra = normalize(spectra)
spectra = aggregate(spectra, meta(spectra)$X, mean, try_keep_txt(mean))
spec_df = as.data.frame(spectra)
data = meta(spectra)
spec.m = as.matrix(spectra) * 100
spectra_percent = as_spectra(spec.m)

#vn_spectra = normalize(spectra_percent)
#vn_spec_df = as.data.frame(vn_spectra)
#data$normalization_magnitude = meta(vn_spectra)$normalization_magnitude

meta(spectra_percent) = data
spec_df = as.data.frame(spectra_percent)



m1_aic = c()
m2_aic = c()

m1_bic = c()
m2_bic = c()

for(i in seq(400, 2400, 1)) {
    x = toString(i)
   
    m1 = lmer(spec_df[, x] ~ age  + (1|scientificName),
                    data = spec_df, REML = T, 
                    lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
    m1_aic = append(m1_aic, AIC(m1))
    m1_bic = append(m1_bic, BIC(m1))
    
    m2 = lmer(spec_df[, x] ~ age + (1 + age|scientificName),
                        data = spec_df, REML = T, 
                        lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
    m2_aic = append(m2_aic, AIC(m2))
    m2_bic = append(m2_bic, BIC(m2))
}

aic_dif = m1_aic - m2_aic

jpeg(filename = '../../lichen figures/aic_comparision_fixed-slope_minus_var-slope_noNorm.jpeg',
     width = 8, height = 6, units = 'in', res = 1200)
par(mfrow = c(1,1))
wv = seq(400, 2400, 1)
plot(wv, aic_dif, type = 'l', xlab = 'Wavelength (nm)', ylab = '∆AIC', 
     main = 'Fixed-slope AIC minus variable-slope AIC')
abline(h=0, col = 'blue', lty = 2)
rect(0, -2, 2500, 2,
        col = rgb(red= 0, green=0, 
                  blue=0, alpha=0.1),
        lty = 0)
dev.off()

plot(wv, m1_bic, col = 'blue', type = 'l', xlab = 'Wavelength (nm)', ylab = 'BIC')
lines(wv, m2_bic, col = 'red')
legend('bottomright',
       c('M1', 'M2'),
       col = c('blue', 'red'), lty = c(1,1))

################################################################################
# Species as a random effect - with vector norm - no normMag
################################################################################
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
spectra = normalize(spectra)
spectra = aggregate(spectra, meta(spectra)$X, mean, try_keep_txt(mean))
spec_df = as.data.frame(spectra)

ranIntercepts = as.data.frame(matrix(nrow = 29))
rownames(ranIntercepts) = sort(unique(spec_df$scientificName))


intercept_variance = c()
resid_variance = c()
psuedoR2 = c()

fixed_intercept_list = c()
fixed_slope_list = c()

intercept_2.5_list = c()
intercept_97.5_list = c()
slope_97.5_list = c()
slope_2.5_list = c()

condR2 = c()
condR2_lower = c()
condR2_upper = c()

marR2 = c()
marR2_lower = c()
marR2_upper = c()


for(i in seq(400, 2400, 1)) {
  x = toString(i)
  lmm = lmer(spec_df[, x] ~ age + (1|scientificName),
             data = spec_df, REML = T, 
             lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
  
  coefs = coef(lmm)
  ranIntercepts = cbind(ranIntercepts, coefs[[1]][,1])
  
  variances = as.data.frame(VarCorr(lmm))
  intercept_variance = append(intercept_variance, variances[1,4])
  resid_variance = append(resid_variance, variances[2,4])
  psuedoR2 = append(psuedoR2, (variances[1,4]/(variances[1,4] + variances[2,4])))
  
  fixed_intercept_list = append(fixed_intercept_list, as.numeric(fixef(lmm)[1]))
  fixed_slope_list = append(fixed_slope_list, as.numeric(fixef(lmm)[2]))
  
  ci = confint(lmm)
  intercept_2.5_list = append(intercept_2.5_list, ci[3])
  intercept_97.5_list = append(intercept_97.5_list, ci[7])
  slope_2.5_list = append(slope_2.5_list, ci[4])
  slope_97.5_list = append(slope_97.5_list, ci[8])
  
  condPart = partR2(lmm, data = spec_df, R2_type = "conditional", nboot = 10)
  condR2 = append(condR2, condPart[3]$R2$estimate)
  condR2_lower = append(condR2_lower, condPart[3]$R2$CI_lower)
  condR2_upper = append(condR2_upper, condPart[3]$R2$CI_upper)
  
  marPart = partR2(lmm, data = spec_df, R2_type = "marginal", nboot = 10)
  marR2 = append(marR2, marPart[3]$R2$estimate)
  marR2_lower = append(marR2_lower, marPart[3]$R2$CI_lower)
  marR2_upper = append(marR2_upper, marPart[3]$R2$CI_upper)
}

stats_list = list()

ranIntercepts = ranIntercepts[,-1]
colnames(ranIntercepts) = seq(400,2400, 1)


stats_list = list.append(stats_list, ranIntercepts) #1
stats_list = list.append(stats_list, intercept_variance) #2
stats_list = list.append(stats_list, resid_variance) #3
stats_list = list.append(stats_list, psuedoR2) #4
stats_list = list.append(stats_list, fixed_intercept_list) #5
stats_list = list.append(stats_list, fixed_slope_list) #6
stats_list = list.append(stats_list, intercept_2.5_list) #7
stats_list = list.append(stats_list, intercept_97.5_list) #8
stats_list = list.append(stats_list, slope_2.5_list) #9
stats_list = list.append(stats_list, slope_97.5_list) #10
stats_list = list.append(stats_list, condR2) #11
stats_list = list.append(stats_list, condR2_lower) #12
stats_list = list.append(stats_list, condR2_upper) #13
stats_list = list.append(stats_list, marR2) # 14
stats_list = list.append(stats_list, marR2_lower) #15
stats_list = list.append(stats_list, marR2_upper) #16


saveRDS(stats_list, 'models/lmms/lmm_60yrs_fixedSlope_vn.rds')


###############################
#Plot
###############################

stats_list = readRDS('models/lmms/lmm_60yrs_fixedSlope_vn.rds')

jpeg(filename = '../../lichen figures/fixedSlope_vn_results.jpeg',
     width = 10, height = 8, units = 'in', res = 1200)
par(mfrow = c(2,3))
#slopes
wv = seq(400, 2400, 1)
plot(wv, stats_list[[6]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (normalized reflectance/year)',
     ylim = c(min(stats_list[[9]]), max(stats_list[[10]])),
     main = 'Fixed slope')
polygon(c(wv, rev(wv)), c(stats_list[[9]], rev(stats_list[[10]])),
        col = 'grey90',
        lty = 0)
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[6]])

#intercepts
wv = seq(400, 2400, 1)
plot(wv, stats_list[[5]],
     type = 'l', 
     col = rgb(0,0,1,1),
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (normalized reflectance)',
     ylim = c(min(stats_list[[1]]), max(stats_list[[1]])),
     main = 'Intercepts')

for (i in 1:nrow(stats_list[[1]])){
  lines(wv, stats_list[[1]][i,], col = 'grey70' )
}
polygon(c(wv, rev(wv)), c(stats_list[[7]], rev(stats_list[[8]])),
        col = rgb(0,0,1,0.2),
        lty = 0)
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[5]], col = rgb(0,0,1,1))

#marginal and conditional variances
wv = seq(400, 2400, 1)
plot(wv, stats_list[[11]],
     type = 'l',
     col = rgb(0,0,1,1),
     xlab = 'Wavelength (nm)',
     ylab = 'Proportion of variance explained',
     ylim = c(0,1),
     main = 'Conditional and Marginal R2')
polygon(c(wv, rev(wv)), c(stats_list[[12]], rev(stats_list[[13]])),
        col = rgb(0,0,1, 0.2),
        lty = 0)
lines(wv, stats_list[[14]], col = rgb(1,0,0,1))
polygon(c(wv, rev(wv)), c(stats_list[[15]], rev(stats_list[[16]])),
        col = rgb(1,0,0, 0.2),
        lty = 0)
legend('topright',
       c('Conditional R2', 'Marginal R2'),
       col = c(rgb(0,0,1,1), rgb(1,0,0,1)), lty = c(1,1))

#variances
wv = seq(400, 2400, 1)
plot(wv, stats_list[[2]],
     ylim = c(min(stats_list[[2]]), max(stats_list[[2]])),
     main = 'Random effects variance',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[3]], col = 'gray')
legend('topright',
       c('Intercept', 'Residual'),
       col = c('blue', 'gray'), lty = c(1,1,1))

#Random Variance explained
wv = seq(400, 2400, 1)
plot(wv, stats_list[[4]],
     ylim = range(stats_list[[4]]),
     main = 'Random effects variance explained by species',
     ylab = 'Proportion of variance explained',
     xlab = 'Wavelength (nm)',
     type = 'l')


dev.off()

################################################################################
# Species as a random effect - with normalization magnitude
################################################################################
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
spectra = aggregate(spectra, meta(spectra)$X, mean, try_keep_txt(mean))
data = meta(spectra)
spec.m = as.matrix(spectra) * 100
spectra_percent = as_spectra(spec.m)
vn_spectra = normalize(spectra)
data$normalization_magnitude = meta(vn_spectra)$normalization_magnitude
meta(spectra_percent) = data
spec_df = as.data.frame(spectra_percent)


ranIntercepts = as.data.frame(matrix(nrow = 29))
rownames(ranIntercepts) = sort(unique(spec_df$scientificName))

intercept_variance = c()
resid_variance = c()
psuedoR2 = c()

fixed_intercept_list = c()
fixed_age_list = c()
fixed_normMag_list = c()

intercept_2.5_list = c()
intercept_97.5_list = c()
age_97.5_list = c()
age_2.5_list = c()
normMag_2.5_list = c()
normMag_97.5_list = c()

condR2 = c()
condR2_lower = c()
condR2_upper = c()

marR2 = c()
marR2_lower = c()
marR2_upper = c()

ageR2 = c()
ageR2_lower = c()
ageR2_upper = c()

normMagR2 = c()
normMagR2_lower = c()
normMagR2_upper = c()


for(i in seq(400, 2400, 1)) {
  x = toString(i)
  lmm = lmer(spec_df[, x] ~ age + normalization_magnitude + (1|scientificName),
             data = spec_df, REML = T, 
             lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
  
  coefs = coef(lmm)
  ranIntercepts = cbind(ranIntercepts, coefs[[1]][,1])
  
  variances = as.data.frame(VarCorr(lmm))
  intercept_variance = append(intercept_variance, variances[1,4])
  resid_variance = append(resid_variance, variances[2,4])
  psuedoR2 = append(psuedoR2, (variances[1,4]/(variances[1,4] + variances[2,4])))
  
  fixed_intercept_list = append(fixed_intercept_list, as.numeric(fixef(lmm)[1]))
  fixed_age_list = append(fixed_age_list, as.numeric(fixef(lmm)[2]))
  fixed_normMag_list = append(fixed_normMag_list, as.numeric(fixef(lmm)[3]))
  
  ci = confint(lmm)
  intercept_2.5_list = append(intercept_2.5_list, ci[3])
  intercept_97.5_list = append(intercept_97.5_list, ci[8])
  age_2.5_list = append(age_2.5_list, ci[4])
  age_97.5_list = append(age_97.5_list, ci[9])
  normMag_2.5_list = append(normMag_2.5_list, ci[5])
  normMag_97.5_list = append(normMag_97.5_list, ci[10])
  
  cond = partR2(lmm, R2_type = 'conditional', nboot = 10)
  mar = partR2(lmm, partvars = c('age', 'normalization_magnitude'),
               R2_type = 'marginal', nboot = 10)
  
  condR2 = append(condR2, cond$R2$estimate)
  condR2_lower = append(condR2_lower, cond$R2$CI_lower)
  condR2_upper = append(condR2_upper, cond$R2$CI_upper)
  marR2 = append(marR2, mar$R2$estimate[1])
  marR2_lower = append(marR2_lower, mar$R2$CI_lower[1])
  marR2_upper = append(marR2_upper, mar$R2$CI_upper[1])
  ageR2 = append(ageR2, mar$R2$estimate[2])
  ageR2_lower = append(ageR2_lower, mar$R2$CI_lower[2])
  ageR2_upper = append(ageR2_upper, mar$R2$CI_upper[2])
  normMagR2 = append(normMagR2, mar$R2$estimate[3])
  normMagR2_lower = append(normMagR2_lower, mar$R2$CI_lower[3])
  normMagR2_upper = append(normMagR2_upper, mar$R2$CI_upper[3])
}

stats_list = list()

ranIntercepts = ranIntercepts[,-1]
colnames(ranIntercepts) = seq(400,2400, 1)


stats_list = list.append(stats_list, ranIntercepts) #1
stats_list = list.append(stats_list, intercept_variance) #2
stats_list = list.append(stats_list, resid_variance) #3
stats_list = list.append(stats_list, psuedoR2) #4
stats_list = list.append(stats_list, fixed_intercept_list) #5
stats_list = list.append(stats_list, fixed_age_list) #6
stats_list = list.append(stats_list, fixed_normMag_list) #7
stats_list = list.append(stats_list, intercept_2.5_list) #8
stats_list = list.append(stats_list, intercept_97.5_list) #9
stats_list = list.append(stats_list, age_2.5_list) #10
stats_list = list.append(stats_list, age_97.5_list) #11
stats_list = list.append(stats_list, normMag_2.5_list) #12
stats_list = list.append(stats_list, normMag_97.5_list) #13
stats_list = list.append(stats_list, condR2) #14
stats_list = list.append(stats_list, condR2_lower) #15
stats_list = list.append(stats_list, condR2_upper) #16
stats_list = list.append(stats_list, marR2) #17
stats_list = list.append(stats_list, marR2_lower) #18
stats_list = list.append(stats_list, marR2_upper) #19
stats_list = list.append(stats_list, ageR2) #20
stats_list = list.append(stats_list, ageR2_lower) #21
stats_list = list.append(stats_list, ageR2_upper) #22
stats_list = list.append(stats_list, normMagR2) #23
stats_list = list.append(stats_list, normMagR2_lower) #24
stats_list = list.append(stats_list, normMagR2_upper) #25

saveRDS(stats_list, 'models/lmms/lmm_60yrs_fixedSlope_normMag.rds')


###############################
#Plot
###############################

stats_list = readRDS('models/lmms/lmm_60yrs_fixedSlope_normMag.rds')

jpeg(filename = '../../lichen figures/fixedSlope_results_withNormMag_R2.jpeg',
     width = 8, height = 8, units = 'in', res = 1200)
par(mfrow = c(3,2))
#slopes - age
wv = seq(400, 2400, 1)
plot(wv, stats_list[[6]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = '% reflectance / year',
     ylim = c(min(stats_list[[10]]), max(stats_list[[11]])),
     main = 'Effect of age')
polygon(c(wv, rev(wv)), c(stats_list[[10]], rev(stats_list[[11]])),
        col = 'grey90',
        lty = 0)
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[6]])

#slopes - normMag
#wv = seq(400, 2400, 1)
#plot(wv, stats_list[[7]],
#     type = 'l', 
#     xlab = 'Wavelength (nm)', 
#     ylab = 'Reflectance / unit normalization magnitude)',
#     ylim = c(min(stats_list[[12]]), max(stats_list[[13]])),
#     main = 'Effect of normalization magnitude')
#polygon(c(wv, rev(wv)), c(stats_list[[12]], rev(stats_list[[13]])),
#        col = 'grey90',
#        lty = 0)
#abline(h = 0, lty = 2, col = 'blue')
#lines(wv, stats_list[[7]])

#intercepts
wv = seq(400, 2400, 1)
plot(wv, stats_list[[5]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = '% reflectance',
     ylim = c(min(stats_list[[1]]), max(stats_list[[1]])),
     main = 'Intercepts')

for (i in 1:nrow(stats_list[[1]])){
  lines(wv, stats_list[[1]][i,], col = 'grey' )
}
polygon(c(wv, rev(wv)), c(stats_list[[8]], rev(stats_list[[9]])),
        col = rgb(0,0,1,0.2),
        lty = 0)
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[5]], col = rgb(0,0,1,1))

#marginal and conditional variances
wv = seq(400, 2400, 1)
plot(wv, stats_list[[14]],
     type = 'l',
     col = rgb(0,0,1,1),
     xlab = 'Wavelength (nm)',
     ylab = 'R-squared',
     ylim = c(0,1),
     main = 'Conditional and Marginal R-squared')
polygon(c(wv, rev(wv)), c(stats_list[[15]], rev(stats_list[[16]])),
        col = rgb(0,0,1, 0.2),
        lty = 0)
lines(wv, stats_list[[17]], col = rgb(1,0,0,1))
polygon(c(wv, rev(wv)), c(stats_list[[18]], rev(stats_list[[19]])),
        col = rgb(1,0,0, 0.2),
        lty = 0)
legend('bottomright',
       c('Conditional R2', 'Marginal R2'),
       col = c(rgb(0,0,1,1), rgb(1,0,0,1)), lty = c(1,1))

#partial R2 - marginal
wv = seq(400, 2400, 1)
plot(wv, stats_list[[20]],
     type = 'l',
     col = rgb(90/256,180/256,17/2562,1),
     xlab = 'Wavelength (nm)',
     ylab = 'R-squared',
     ylim = c(0,1),
     main = 'Partial R-squared')
polygon(c(wv, rev(wv)), c(stats_list[[21]], rev(stats_list[[22]])),
        col = rgb(90/256,180/256,172/256, 0.2),
        lty = 0)
lines(wv, stats_list[[23]], col = rgb(216/256,179/256,101/256,1))
polygon(c(wv, rev(wv)), c(stats_list[[24]], rev(stats_list[[25]])),
        col = rgb(216/256,179/256,101/256, 0.2),
        lty = 0)
legend('right',
       c('Age R2', 'Normalization magnitude R2'),
       col = c(rgb(90/256,180/256,172/256,1), rgb(216/256,179/256,101/256,1)),
       lty = c(1,1))


#Random variances
wv = seq(400, 2400, 1)
plot(wv, stats_list[[2]],
     ylim = c(min(stats_list[[2]]), max(stats_list[[2]])),
     main = 'Random effects variance',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[3]], col = 'gray')
legend('topright',
       c('Intercept', 'Residual'),
       col = c('blue', 'gray'), lty = c(1,1,1))

#Random Variance explained
wv = seq(400, 2400, 1)
plot(wv, stats_list[[4]],
     ylim = range(stats_list[[4]]),
     main = 'Random effects variance explained by species',
     ylab = 'Proportion of variance explained',
     xlab = 'Wavelength (nm)',
     type = 'l')

dev.off()

################################################################################
# Species as a random effect - no vector normalization
################################################################################
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
spectra = aggregate(spectra, meta(spectra)$X, mean, try_keep_txt(mean))
data = meta(spectra)
spec.m = as.matrix(spectra) * 100
spectra_percent = as_spectra(spec.m)
meta(spectra_percent) = data
spec_df = as.data.frame(spectra_percent)

ranIntercepts = as.data.frame(matrix(nrow = 29))
rownames(ranIntercepts) = sort(unique(spec_df$scientificName))


intercept_variance = c()
resid_variance = c()
psuedoR2 = c()

fixed_intercept_list = c()
fixed_slope_list = c()

intercept_2.5_list = c()
intercept_97.5_list = c()
slope_97.5_list = c()
slope_2.5_list = c()

condR2 = c()
condR2_lower = c()
condR2_upper = c()

marR2 = c()
marR2_lower = c()
marR2_upper = c()

for(i in seq(400, 2400, 1)) {
    x = toString(i)
    lmm = lmer(spec_df[, x] ~ age + (1|scientificName),
                  data = spec_df, REML = T, 
                  lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
    
    coefs = coef(lmm)
    ranIntercepts = cbind(ranIntercepts, coefs[[1]][,1])

    variances = as.data.frame(VarCorr(lmm))
    intercept_variance = append(intercept_variance, variances[1,4])
    resid_variance = append(resid_variance, variances[2,4])
    psuedoR2 = append(psuedoR2, (variances[1,4]/(variances[1,4] + variances[2,4])))
    
    fixed_intercept_list = append(fixed_intercept_list, as.numeric(fixef(lmm)[1]))
    fixed_slope_list = append(fixed_slope_list, as.numeric(fixef(lmm)[2]))
    
    ci = confint(lmm)
    intercept_2.5_list = append(intercept_2.5_list, ci[3])
    intercept_97.5_list = append(intercept_97.5_list, ci[7])
    slope_2.5_list = append(slope_2.5_list, ci[4])
    slope_97.5_list = append(slope_97.5_list, ci[8])
    
    cond = partR2(lmm, R2_type = 'conditional', nboot = 10)
    mar = partR2(lmm, R2_type = 'marginal', nboot = 10)
    
    condR2 = append(condR2, cond$R2$estimate)
    condR2_lower = append(condR2_lower, cond$R2$CI_lower)
    condR2_upper = append(condR2_upper, cond$R2$CI_upper)
    marR2 = append(marR2, mar$R2$estimate[1])
    marR2_lower = append(marR2_lower, mar$R2$CI_lower)
    marR2_upper = append(marR2_upper, mar$R2$CI_upper)
}

stats_list = list()

ranIntercepts = ranIntercepts[,-1]
colnames(ranIntercepts) = seq(400,2400, 1)


stats_list = list.append(stats_list, ranIntercepts) #1
stats_list = list.append(stats_list, intercept_variance) #2
stats_list = list.append(stats_list, resid_variance) #3
stats_list = list.append(stats_list, psuedoR2) #4
stats_list = list.append(stats_list, fixed_intercept_list) #5
stats_list = list.append(stats_list, fixed_slope_list) #6
stats_list = list.append(stats_list, intercept_2.5_list) #7
stats_list = list.append(stats_list, intercept_97.5_list) #8
stats_list = list.append(stats_list, slope_2.5_list) #9
stats_list = list.append(stats_list, slope_97.5_list) #10
stats_list = list.append(stats_list, condR2) #11
stats_list = list.append(stats_list, condR2_lower) #12
stats_list = list.append(stats_list, condR2_upper) #13
stats_list = list.append(stats_list, marR2) #14
stats_list = list.append(stats_list, marR2_lower) #15
stats_list = list.append(stats_list, marR2_upper) #16

saveRDS(stats_list, 'models/lmms/lmm_60yrs_fixedSlope.rds')


###############################
#Plot
###############################

stats_list = readRDS('models/lmms/lmm_60yrs_fixedSlope.rds')

jpeg(filename = '../../lichen figures/fixedSlope_results_noNorm.jpeg',
 width = 8, height = 12, units = 'in', res = 1200)
par(mfrow = c(3,2))
#slopes
wv = seq(400, 2400, 1)
plot(wv, stats_list[[6]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (% reflectance/year)',
     ylim = c(min(stats_list[[9]]), max(stats_list[[10]])),
     main = 'Fixed slope')
polygon(c(wv, rev(wv)), c(stats_list[[9]], rev(stats_list[[10]])),
        col = 'grey90',
        lty = 0)
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[6]])

#intercepts
wv = seq(400, 2400, 1)
plot(wv, stats_list[[5]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (% reflectance)',
     ylim = c(min(stats_list[[1]]), max(stats_list[[1]])),
     main = 'Intercepts')
for (i in 1:nrow(stats_list[[1]])){
  lines(wv, stats_list[[1]][i,], col = 'grey' )
}
polygon(c(wv, rev(wv)), c(stats_list[[7]], rev(stats_list[[8]])),
        col = rgb(0,0,1,0.2),
        lty = 0)
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[5]], col = rgb(0,0,0,1))

#marginal and conditional variances
wv = seq(400, 2400, 1)
plot(wv, stats_list[[11]],
     type = 'l',
     col = rgb(0,0,1,1),
     xlab = 'Wavelength (nm)',
     ylab = 'Proportion of variance explained',
     ylim = c(0,1),
     main = 'R-squared')
polygon(c(wv, rev(wv)), c(stats_list[[12]], rev(stats_list[[13]])),
        col = rgb(0,0,1, 0.2),
        lty = 0)
lines(wv, stats_list[[14]], col = rgb(1,0,0,1))
polygon(c(wv, rev(wv)), c(stats_list[[15]], rev(stats_list[[16]])),
        col = rgb(1,0,0, 0.2),
        lty = 0)
legend('topright',
       c('Conditional R2', 'Marginal R2'),
       col = c(rgb(0,0,1,1), rgb(1,0,0,1)), lty = c(1,1))


#variances
wv = seq(400, 2400, 1)
plot(wv, stats_list[[2]],
     ylim = c(min(stats_list[[2]]), max(stats_list[[2]])),
     main = 'Random effects variance',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[3]], col = 'gray')
legend('topright',
       c('Intercept', 'Residual'),
       col = c('blue', 'gray'), lty = c(1,1,1))

#Random Variance explained
wv = seq(400, 2400, 1)
plot(wv, stats_list[[4]],
     ylim = range(stats_list[[4]]),
     main = 'Random effects variance explained by species',
     ylab = 'Proportion of variance explained',
     xlab = 'Wavelength (nm)',
     type = 'l')

dev.off()

################################################################################
# Species as a random effect - with vector normalization
################################################################################
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = normalize(spectra[meta(spectra)$age <= 60, ])
spectra = aggregate(spectra, meta(spectra)$X, mean, try_keep_txt(mean))
spec_df = as.data.frame(spectra)

ranIntercepts = as.data.frame(matrix(nrow = 29))
rownames(ranIntercepts) = sort(unique(spec_df$scientificName))


intercept_variance = c()
resid_variance = c()
psuedoR2 = c()

fixed_intercept_list = c()
fixed_age_list = c()
fixed_normMag_list = c()

intercept_2.5_list = c()
intercept_97.5_list = c()
age_97.5_list = c()
age_2.5_list = c()
normMag_2.5_list = c()
normMag_97.5_list = c()


for(i in seq(400, 2400, 1)) {
  x = toString(i)
  lmm = lmer(spec_df[, '700'] ~ age + normalization_magnitude + (1|scientificName),
             data = spec_df, REML = T, 
             lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
  
  coefs = coef(lmm)
  ranIntercepts = cbind(ranIntercepts, coefs[[1]][,1])
  
  variances = as.data.frame(VarCorr(lmm))
  intercept_variance = append(intercept_variance, variances[1,4])
  resid_variance = append(resid_variance, variances[2,4])
  psuedoR2 = append(psuedoR2, (variances[1,4]/(variances[1,4] + variances[2,4])))
  
  fixed_intercept_list = append(fixed_intercept_list, as.numeric(fixef(lmm)[1]))
  fixed_age_list = append(fixed_age_list, as.numeric(fixef(lmm)[2]))
  fixed_normMag_list = append(fixed_normMag_list, as.numeric(fixef(lmm)[3]))
  
  ci = confint(lmm)
  intercept_2.5_list = append(intercept_2.5_list, ci[3])
  intercept_97.5_list = append(intercept_97.5_list, ci[8])
  age_2.5_list = append(age_2.5_list, ci[4])
  age_97.5_list = append(age_97.5_list, ci[9])
  normMag_2.5_list = append(normMag_2.5_list, ci[5])
  normMag_97.5_list = append(normMag_97.5_list, ci[10])
}

stats_list = list()

ranIntercepts = ranIntercepts[,-1]
colnames(ranIntercepts) = seq(400,2400, 1)


stats_list = list.append(stats_list, ranIntercepts) #1
stats_list = list.append(stats_list, intercept_variance) #2
stats_list = list.append(stats_list, resid_variance) #3
stats_list = list.append(stats_list, psuedoR2) #4
stats_list = list.append(stats_list, fixed_intercept_list) #5
stats_list = list.append(stats_list, fixed_age_list) #6
stats_list = list.append(stats_list, fixed_normMag_list) #7
stats_list = list.append(stats_list, intercept_2.5_list) #8
stats_list = list.append(stats_list, intercept_97.5_list) #9
stats_list = list.append(stats_list, age_2.5_list) #10
stats_list = list.append(stats_list, age_97.5_list) #11
stats_list = list.append(stats_list, normMag_2.5_list) #12
stats_list = list.append(stats_list, normMag_97.5_list) #13

saveRDS(stats_list, 'models/lmms/lmm_60yrs_fixedSlope_vn.rds')


###############################
#Plot
###############################

stats_list = readRDS('models/lmms/lmm_60yrs_fixedSlope_vn.rds')

jpeg(filename = '../../lichen figures/fixedSlope_results_vn.jpeg',
     width = 12, height = 8, units = 'in', res = 1200)
par(mfrow = c(2,3))
#slopes - age
wv = seq(400, 2400, 1)
plot(wv, stats_list[[6]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Vector normalized reflectance / year',
     ylim = c(min(stats_list[[10]]), max(stats_list[[11]])),
     main = 'Effect of age')
polygon(c(wv, rev(wv)), c(stats_list[[10]], rev(stats_list[[11]])),
        col = 'grey90',
        lty = 0)
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[6]])

#slopes - normMag
wv = seq(400, 2400, 1)
plot(wv, stats_list[[7]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Reflectance / unit normalization magnitude)',
     ylim = c(min(stats_list[[12]]), max(stats_list[[13]])),
     main = 'Effect of normalization magnitude')
polygon(c(wv, rev(wv)), c(stats_list[[12]], rev(stats_list[[13]])),
        col = 'grey90',
        lty = 0)
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[7]])

#intercepts
wv = seq(400, 2400, 1)
plot(wv, stats_list[[5]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Vector normalized reflectance',
     ylim = c(min(stats_list[[1]]), max(stats_list[[1]])),
     main = 'Intercepts')
polygon(c(wv, rev(wv)), c(stats_list[[8]], rev(stats_list[[9]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[1]])){
  lines(wv, stats_list[[1]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[5]])

#variances
wv = seq(400, 2400, 1)
plot(wv, stats_list[[2]],
     ylim = c(min(stats_list[[2]]), max(stats_list[[2]])),
     main = 'Random effects variance',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[3]], col = 'gray')
legend('topright',
       c('Intercept', 'Residual'),
       col = c('blue', 'gray'), lty = c(1,1,1))

#Random Variance explained
wv = seq(400, 2400, 1)
plot(wv, stats_list[[4]],
     ylim = range(stats_list[[4]]),
     main = 'Random effects variance explained by species',
     ylab = 'Proportion of variance explained',
     xlab = 'Wavelength (nm)',
     type = 'l')

dev.off()

################################################################################
#Brightness as a function of age
################################################################################
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
spectra = normalize(spectra)
spectra = aggregate(spectra, meta(spectra)$X, mean, try_keep_txt(mean))
spec_df = as.data.frame(spectra)

m = lm(normalization_magnitude ~ age, data = spec_df)
lmm = lmer(normalization_magnitude ~ age + (1|scientificName), data = spec_df)
lmm2 = lmer(normalization_magnitude ~ age + (1 + age|scientificName),
            data = spec_df, REML = T, 
            lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5,
                        optCtrl = list(maxfun = 1e5)))

AIC(m)
AIC(lmm)
AIC(lmm2)

partR2(lmm, data = spec_df, R2_type = "marginal", nboot = 10)
partR2(lmm, data = spec_df, R2_type = "conditional", nboot = 10)

#plot
coefs = coef(lmm)

jpeg(filename = '../../lichen figures/normMag-age_allspecies.jpeg',
     width = 9, height = 8, units = 'in', res = 1200)
plot(spec_df$age, spec_df$normalization_magnitude, ylab = 'Normalization Magnitude',
     xlab = ' Age')
for (i in 1:29) {
  abline(coefs$scientificName$`(Intercept)`[i], coefs$scientificName$age[i],
         col = 'gray50')
}
abline(fixef(lmm)[1], fixef(lmm)[2], col = 'blue', lwd = 2)
dev.off()

library(ggeffects)
pred <- ggpredict(lmm, terms = c("age"))

jpeg(filename = '../../lichen figures/normMag-age.jpeg',
     width = 8, height = 3, units = 'in', res = 1200)
(ggplot(pred) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = spec_df,                      
               aes(x = age, y = normalization_magnitude)) + #, colour = scientificName)) + 
    labs(x = "Age", y = "Normalization magnitude", 
         title = "") + 
    theme_minimal()
)

dev.off()

ggpredict(lmm, terms = c("age", "scientificName"), type = "re") %>% 
  plot() +
  labs(x = "Age", y = "Normalization Magnitude", title = "") + 
  scale_fill_manual(values = rep('', 29)) +
  theme_minimal()


library(sjPlot)
# Visualise random effects 
(re.effects <- plot_model(lmm, type = "re", show.values = TRUE))


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
stats_list = readRDS('models/lmms/lmm_hier2_60yrs.rds')


################################################################################
#Family
par(mfrow = c(2,1))
##slopes
wv = seq(400, 2400, 10)
plot(wv, stats_list[[15]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (% reflectance/year)',
     ylim = c(min(stats_list[[4]]), max(stats_list[[19]])),
     main = 'Slopes - Family')
polygon(c(wv, rev(wv)), c(stats_list[[18]], rev(stats_list[[19]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[4]])){
  lines(wv, stats_list[[4]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[15]])

##intercepts
wv = seq(400, 2400, 10)
plot(wv, stats_list[[14]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (% reflectance)',
     ylim = c(min(stats_list[[1]]), max(stats_list[[17]])),
     main = 'Intercepts - Family')
polygon(c(wv, rev(wv)), c(stats_list[[16]], rev(stats_list[[17]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[1]])){
  lines(wv, stats_list[[1]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[14]])

################################################################################
#Order
par(mfrow = c(2,1))
##slopes
wv = seq(400, 2400, 10)
plot(wv, stats_list[[15]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (% reflectance/year)',
     ylim = c(min(stats_list[[18]]), max(stats_list[[19]])),
     main = 'Slopes - Order')
polygon(c(wv, rev(wv)), c(stats_list[[18]], rev(stats_list[[19]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[5]])){
  lines(wv, stats_list[[5]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[15]])

##intercepts
wv = seq(400, 2400, 10)
plot(wv, stats_list[[14]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (% reflectance)',
     ylim = c(min(stats_list[[2]]), max(stats_list[[17]])),
     main = 'Intercepts - Order')
polygon(c(wv, rev(wv)), c(stats_list[[16]], rev(stats_list[[17]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[2]])){
  lines(wv, stats_list[[2]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[14]])

################################################################################
#Class
par(mfrow = c(2,1))
##slopes
wv = seq(400, 2400, 10)
plot(wv, stats_list[[15]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (% reflectance/year)',
     ylim = c(min(stats_list[[6]]), max(stats_list[[6]])),
     main = 'Slopes - Class')
polygon(c(wv, rev(wv)), c(stats_list[[18]], rev(stats_list[[19]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[6]])){
  lines(wv, stats_list[[6]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[15]])

##intercepts
wv = seq(400, 2400, 10)
plot(wv, stats_list[[14]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Intercept (% reflectance)',
     ylim = c(min(stats_list[[16]]), max(stats_list[[3]])),
     main = 'Intercepts - Class')
polygon(c(wv, rev(wv)), c(stats_list[[16]], rev(stats_list[[17]])),
        col = 'grey90',
        lty = 0)
for (i in 1:nrow(stats_list[[3]])){
  lines(wv, stats_list[[3]][i,], col = 'grey' )
}
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[14]])

################################################################################
#variances
par(mfrow = c(3,2))

##Family
wv = seq(400, 2400, 10)
plot(wv, stats_list[[7]],
     ylim = c(0, max(stats_list[[7]])),
     main = 'Variance - Family',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[10]], col = 'red')
lines(wv, stats_list[[13]], col = 'gray')

plot(wv, stats_list[[10]], main = 'Slope variance - Family', ylab = 'variance',
     xlab = 'Wavelength (nm)', type = 'l')

##Order
wv = seq(400, 2400, 10)
plot(wv, stats_list[[8]],
     ylim = c(0, max(stats_list[[13]])),
     main = 'Variance - Order',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[11]], col = 'red')
lines(wv, stats_list[[13]], col = 'gray')

plot(wv, stats_list[[11]], main = 'Slope variance - Order', ylab = 'variance',
     xlab = 'Wavelength (nm)', type = 'l')

##Class
wv = seq(400, 2400, 10)
plot(wv, stats_list[[9]],
     ylim = c(min(stats_list[[9]]), max(stats_list[[9]])),
     main = 'Variance - Class',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = 'blue',
     type = 'l')
lines(wv, stats_list[[12]], col = 'red')
lines(wv, stats_list[[13]], col = 'gray')

plot(wv, stats_list[[12]], main = 'Slope variance - Class', ylab = 'variance',
     xlab = 'Wavelength (nm)', type = 'l')

##Total slope variance
par(mfrow = c(1,1))
wv = seq(400, 2400, 10)
plot(wv, stats_list[[10]],
     ylim = c(min(stats_list[[10]], stats_list[[11]], stats_list[[12]]),
              max(stats_list[[10]], stats_list[[11]], stats_list[[12]])),
     main = 'Variance - Slopes',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = '#a6cee3',
     type = 'l')
lines(wv, stats_list[[11]], col = '#1f78b4')
lines(wv, stats_list[[12]], col = '#b2df8a')
legend('topright',
       c('Family', 'Order', 'Class'),
       col = c('#a6cee3', '#1f78b4', '#b2df8a'), lty = c(1,1,1,1))

##Total intercept variance
wv = seq(400, 2400, 10)
plot(wv, stats_list[[7]],
     ylim = c(min(stats_list[[7]], stats_list[[8]], stats_list[[9]]),
              max(stats_list[[7]], stats_list[[8]], stats_list[[9]])),
     main = 'Variance - Intercepts',
     ylab = 'Variance',
     xlab = 'Wavelength (nm)',
     col = '#a6cee3',
     type = 'l')
lines(wv, stats_list[[8]], col = '#1f78b4')
lines(wv, stats_list[[9]], col = '#b2df8a')
legend('bottomright',
       c('Family', 'Order', 'Class'),
       col = c('#a6cee3', '#1f78b4', '#b2df8a'), lty = c(1,1,1,1))

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

