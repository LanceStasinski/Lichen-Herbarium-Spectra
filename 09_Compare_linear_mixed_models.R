################################################################################
#Compare models - intercept only, fixed slope, variable slope
################################################################################
library(spectrolab)
library(lme4)
library(nlme)
library(rlist)
library(optimx)
library(partR2)

#setup spectra
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
#flavosToRemove = c('Flavoparmelia_baltimorensis', 'Flavoparmelia_euplecta', 
#                   'Flavoparmelia_haysomii', 'Flavoparmelia_rutidota', 
#                   'Flavoparmelia_soredians', 'Flavopunctelia_flaventior',
#                   'Flavopunctelia_praesignis', 'Flavopunctelia_soredica')

#for (i in seq(1, length(flavosToRemove), 1)) {
#  spectra = spectra[!meta(spectra)$scientificName == flavosToRemove[i]]
#}

spectra = spectra[!meta(spectra)$Morphology == 'fruticose',]
spectra = spectra[!meta(spectra)$Morphology == 'crustose',]

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


#create empty vectors
m1_aic = c()
m2_aic = c()

m1_bic = c()
m2_bic = c()

m3_aic = c()
m3_bic = c()

#run models and append AIC and BIC metrics
for(i in seq(400, 2400, 1)) {
  x = toString(i)
  
  m1 = lmer(spec_df[, x] ~ age  + (1|scientificName),
            data = spec_df, REML = T, 
            lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
  m1_aic = append(m1_aic, AIC(m1))
  m1_bic = append(m1_bic, BIC(m1))
  
  m2 = lmer(spec_df[, x] ~ 1 + (1|scientificName),
            data = spec_df, REML = T, 
            lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
  m2_aic = append(m2_aic, AIC(m2))
  m2_bic = append(m2_bic, BIC(m2))
  
  m3 = lmer(spec_df[, x] ~ 1 + (1 + age|scientificName),
            data = spec_df, REML = T, 
            lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5, optCtrl = list(maxfun = 1e5)))
  m3_aic = append(m3_aic, AIC(m3))
  m3_bic = append(m3_bic, BIC(m3))
  
}

#calculate Delta AIC (fixed slope AIC minus AIC of other models)
aic_dif1 = m1_aic - m2_aic
aic_dif2 = m1_aic - m3_aic

#plot delta AIC
jpeg(filename = '../../lichen figures/AIC/species_foliose_delta_AIC.jpeg',
     width = 8, height = 6, units = 'in', res = 1200)
par(mfrow = c(1,1))
wv = seq(400, 2400, 1)
plot(wv, aic_dif1, type = 'l', xlab = 'Wavelength (nm)', ylab = '∆AIC', 
     main = 'Species as a random effect - Foliose', col = 'blue', ylim = c(-5,15))
abline(h=0, col = 'black', lty = 2)
rect(0, -2, 2500, 2,
     col = rgb(red= 0, green=0, 
               blue=0, alpha=0.1),
     lty = 0)
lines(wv, aic_dif2, col = 'red')
legend('bottomright',
       c('∆AIC Intercept-only', '∆AIC Variable-slope'),
       col = c('blue', 'red'), lty = c(1,1))
dev.off()
