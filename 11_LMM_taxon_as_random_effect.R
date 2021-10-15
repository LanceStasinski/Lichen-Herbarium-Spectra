################################################################################
# Species (or other taxon) as a random effect - no vector normalization
################################################################################
library(spectrolab)
library(lme4)
library(nlme)
library(rlist)
library(optimx)
library(partR2)

#setup spectral data
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
spectra = aggregate(spectra, meta(spectra)$X, mean, try_keep_txt(mean))
data = meta(spectra)
spec.m = as.matrix(spectra) * 100
spectra_percent = as_spectra(spec.m)
meta(spectra_percent) = data
spec_df = as.data.frame(spectra_percent)

#setup empty vectors
ranIntercepts = as.data.frame(matrix(nrow = 29))
rownames(ranIntercepts) = sort(unique(spec_df$Family))

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

#run model for each wavelength and append metrics
for(i in seq(400, 2400, 1)) {
  x = toString(i)
  lmm = lmer(spec_df[, x] ~ age + (1|Family),
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

#combine metrics into a single list
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

saveRDS(stats_list, 'models/lmms/lmm_60yrs_fixedSlope_family.rds')


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
     ylab = 'r-squared',
     ylim = c(0,1),
     main = 'Conditional and Marginal r-squared')
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