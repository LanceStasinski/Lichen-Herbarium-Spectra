#Linear Mixed Models
################################################################################
#setup
################################################################################

library(spectrolab)

################################################################################
#comparison of functions - testing
################################################################################

spectra = readRDS('spectra/lichen_spectra.rds')
spec_df = as.data.frame(spectra)

#Standardize explanatory variable (age)
spec_df$age = scale(spec_df$age, center = T, scale = T)

#basic linear model
basic.lm <- lm(`700` ~ age, data = spec_df)
summary(basic.lm)

#plot
library(ggplot2)

(prelim_plot <- ggplot(spec_df, aes(x = age, y = `700`)) +
    geom_point() +
    geom_smooth(method = "lm"))

plot(basic.lm, which = 1) #Residuals vs fitted
plot(basic.lm, which = 2) #QQ plot

#lme
library(Rcpp)
library(lme4)
mixed.lmer <- lmer(`700` ~ age + (1|scientificName), data = spec_df)
summary(mixed.lmer)



