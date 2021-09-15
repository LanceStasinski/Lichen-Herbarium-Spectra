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
    lmm_varSlope = lmer(spec_df[, '2050'] ~ age + (age|scientificName),
                        data = spec_df, REML = T, 
                        lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    lmm_fixedSlope = lmer(spec_df[, x] ~ age + (1|Class/Order/Family), data = spec_df, REML = T, 
                         control = lmerControl(
                           optimizer ='Nelder_Mead'))
    
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
################################################################################
#Bayesian approach
################################################################################
library(blme)
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]

spec_df = as.data.frame(spectra)
spec_df$age = scale(spec_df$age, center = TRUE, scale = TRUE)

lmm = blmer(spec_df[, '700'] ~ age + (age|scientificName), data = spec_df, REML = T, 
            control = lmerControl(
                optimizer ='optimx', optCtrl=list(method='nlminb')))









