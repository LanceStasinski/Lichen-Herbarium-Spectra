################################################################################
#Brightness as a function of age
################################################################################

library(spectrolab)
library(lme4)
library(nlme)
library(rlist)
library(optimx)
library(partR2)


#setup spectra
spectra = readRDS('spectra/lichen_spectra.rds')
#400:700
#701:1100
#1101:2400
spectra = spectra[,700:1900 ]
spectra = spectra[meta(spectra)$age <= 60, ]
spectra = normalize(spectra)
spectra = aggregate(spectra, meta(spectra)$X, mean, try_keep_txt(mean))
spec_df = as.data.frame(spectra)

m = lm(normalization_magnitude ~ age, data = spec_df)
intercept = lmer(normalization_magnitude ~ 1 + (1|scientificName), data = spec_df)
fixed_slope = lmer(normalization_magnitude ~ age + (1|scientificName), data = spec_df)
var_slope = lmer(normalization_magnitude ~ age + (1 + age|scientificName),
                 data = spec_df, REML = T, 
                 lmerControl(optimizer ='bobyqa', boundary.tol = 1e-5,
                             optCtrl = list(maxfun = 1e5)))


AIC(m)
AIC(intercept)
AIC(fixed_slope)
AIC(var_slope)


partR2(lmm, data = spec_df, R2_type = "marginal", nboot = 10)
partR2(lmm, data = spec_df, R2_type = "conditional", nboot = 10)

#####################################
#plot
#####################################

#Colors
library(RColorBrewer)
col_pal = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, col_pal$maxcolors, rownames(col_pal)))
cols = sample(col_vector, 29)
spec_df$color = 'black'
species = sort(unique(spec_df$scientificName))

for (i in 1:length(species)) {
  spec_df$color[spec_df$scientificName == species[i]] = cols[i]
}


coefs = coef(fixed_slope)

jpeg(filename = '../../lichen figures/normMag-age_700-1900.jpeg',
     width = 9, height = 8, units = 'in', res = 1200)
plot(spec_df$age, spec_df$normalization_magnitude, ylab = 'Normalization Magnitude',
     xlab = ' Age', col = spec_df$color, pch = 16)
for (i in 1:29) {
  abline(coefs$scientificName$`(Intercept)`[i], coefs$scientificName$age[i],
         col = cols[i])
}
abline(fixef(fixed_slope)[1], fixef(fixed_slope)[2], col = 'black', lwd = 3)
dev.off()
