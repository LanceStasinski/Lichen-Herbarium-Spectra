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

spec_df = as.data.frame(spectra)
spec_df$age = scale(spec_df$age, center = TRUE, scale = TRUE)


varSlope_aic = c()
fixedSlope_aic = c()

for(i in seq(400, 2400, 10)) {
    x = toString(i)
    lmm_varSlope = lmer(spec_df[, x] ~ age + (age|Class), data = spec_df, REML = T, 
                        lmerControl(
                            optimizer ='optimx', optCtrl=list(method='nlminb')))
    lmm_fixedSlope = lmer(spec_df[, x] ~ age + (1|Class), data = spec_df, REML = T, 
                         control = lmerControl(
                             optimizer ='optimx', optCtrl=list(method='nlminb')))
    
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
spec_df = as.data.frame(spectra)
#spec_df$age = scale(spec_df$age, center = TRUE, scale = TRUE)

ranEffects = as.data.frame(matrix(nrow = 16))
rownames(ranEffects) = sort(unique(spec_df$Order))
int_r2_list = c()
slope_r2_list = c()
age_effect_list = c()
age_97.5_list = c()
age_2.5_list = c()


for(i in seq(400, 2400, 1)) {
    x = toString(i)
    lmm = lmer(spec_df[, x] ~ age + (age|Family), data = spec_df, REML = T, 
               control = lmerControl(
                   optimizer ='optimx', optCtrl=list(method='nlminb')))
    
    ranEf = ranef(lmm)
    ranEffects = cbind(ranEffects, as.data.frame(ranEf[[1]][2]))
    d = as.data.frame(VarCorr(lmm))
    intercept_var = d[1,4]
    slope_var = d[2,4]
    residual_var = d[4,4]
    total_var = sum(intercept_var, slope_var, residual_var)
    
    int_r2 = intercept_var/total_var
    int_r2_list = append(int_r2_list, int_r2)
    slope_r2 = slope_var/total_var
    slope_r2_list = append(slope_r2_list, slope_r2)
    
    age_effect_list = append(age_effect_list, as.numeric(fixef(lmm)[2]))
    
    
    ci = confint(lmm, method = 'Wald')
    age_2.5_list = append(age_2.5_list, ci[6])
    age_97.5_list = append(age_97.5_list, ci[12])
}

stats_list = list()
stats_list = list.append(stats_list, int_r2_list)
stats_list = list.append(stats_list, slope_r2_list)
stats_list = list.append(stats_list, age_effect_list)
stats_list = list.append(stats_list, age_97.5_list)
stats_list = list.append(stats_list, age_2.5_list)
ranEffects = ranEffects[,-1]
colnames(ranEffects) = seq(400, 2400, 1)
stats_list = list.append(stats_list, ranEffects)

saveRDS(stats_list, 'models/lmm_60yrs_family.rds')

stats_list = readRDS('models/lmm_1_vn.rds')

par(mfrow = c(1,1))
wv = seq(400, 2400, 1)
plot(wv, stats_list[[3]],
     type = 'l', 
     xlab = 'Wavelength (nm)', 
     ylab = 'Effect of age (age/scaled(wavelength))',
     ylim = c(min(stats_list[[5]]), max(stats_list[[4]])))
polygon(c(wv, rev(wv)), c(stats_list[[4]], rev(stats_list[[5]])),
        col = 'grey80',
        lty = 0)
abline(h = 0, lty = 2, col = 'blue')
lines(wv, stats_list[[3]])


plot(wv, stats_list[[1]],
     type = 'l',
     ylab = 'Percent variance explained',
     xlab = 'Wavelength (nm)',
     col = 'red',
     ylim = c(0, 1))
lines(wv, stats_list[[2]], col = 'blue')
legend('topright',
       c('Intercept', 'Slope'),
       col = c('red', 'blue'), lty = c(1,1))


