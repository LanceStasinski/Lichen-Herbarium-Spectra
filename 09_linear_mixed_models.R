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

lmm_varslope = lmer(`700` ~ age + (age|scientificName), data = spec_df, REML = FALSE, 
                    +                     control = lmerControl(optimizer ="Nelder_Mead"))

################################################################################
#comparison of functions - testing
################################################################################

spectra = readRDS('spectra/lichen_spectra.rds')
#4 scans per individual cannot be treated as independent in a linear model, so 
#I'm taking the mean per individual
spectra = aggregate(spectra, by = meta(spectra)$X, mean, try_keep_txt(mean))
spec_df = as.data.frame(spectra)

#Standardize explanatory variable (age)
spec_df$age = scale(spec_df$age, center = T, scale = T)

#basic linear model
basic.lm <- lm(`700` ~ age, data = spec_df)
summary(basic.lm)

#plot

(prelim_plot <- ggplot(spec_df, aes(x = age, y = `700`)) +
    geom_point() +
    geom_smooth(method = "lm"))

plot(basic.lm, which = 1) #Residuals vs fitted
plot(basic.lm, which = 2) #QQ plot

#lme
mixed.lmer <- lmer(`700` ~ age + (1|scientificName), data = spec_df)
summary(mixed.lmer)
#Species explains 62% of the variance in reflectance that is left over after the variance explained by age.  
#Absolute value of slope estimates is larger than the error 
#which indicates the effect is distinguishable from 0.
plot(mixed.lmer)
qqnorm(resid(mixed.lmer))
qqline(resid(mixed.lmer))


#variable intercept, fixed slope - nested
mixed.lmer2 = lmer(`700` ~ age + (1|Class/Order/Family), data = spec_df)
summary(mixed.lmer2)

#random slope - random intercept
mixed.ranslope = lmer(`700` ~ age + (1 + age|scientificName), data = spec_df)
summary(mixed.ranslope)
anova(mixed.ranslope)
#Plot more elegantly

# Extract the prediction data frame
pred.mm <- ggpredict(mixed.ranslope, terms = c("age"))  # this gives overall predictions for the model

# Plot the predictions 

(ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = spec_df,                      # adding the raw data (scaled values)
               aes(x = age, y = `700`, colour = scientificName)) + 
    labs(x = "Age (scaled)", y = "Reflectance", 
         title = "Reflectance vs age") + 
    theme_minimal()
)



# Visualise random effects 
(re.effects <- plot_model(mixed.ranslope, type = "re", show.values = TRUE))



################################################################################

spectra = readRDS('spectra/lichen_spectra.rds')
#4 scans per individual cannot be treated as independent in a linear model, so 
#I'm taking the mean per individual
spectra = aggregate(spectra, by = meta(spectra)$X, mean, try_keep_txt(mean))
spec_df = as.data.frame(spectra)

uniqueNames = unique(spec_df$scientificName)
nameList = c()
for (i in 1:length(uniqueNames)) {
    speciesSpec = spec_df[spec_df$scientificName == uniqueNames[i],]
    if (nrow(speciesSpec) > 4) {
        nameList = append(nameList, uniqueNames[i])
    }
}
spec_df = spec_df[spec_df$scientificName %in% nameList, ]

linear = lm(`700` ~ age, data = spec_df)
glm = lm(`700` ~ age + scientificName, data = spec_df)
lmm = lme(`700`~ age, data = spec_df, random = ~1|scientificName)
lmm_varslope = lmer(`700` ~ age + (age|scientificName), data = spec_df)
lmm_varslope_sqrt = lmer(`700` ~ sqrt(age) + (sqrt(age)|scientificName), data = spec_df)
lmm_varslope_log10 = lmer(`700` ~ log10(age) + (log10(age)|scientificName), data = spec_df)
lmm_varslope_log = lmer(`700` ~ log(age) + (log(age)|scientificName), data = spec_df)

lmm_varslope_log10_2 = lmer(log10(`700`) ~ log10(age) + (log10(age)|scientificName), data = spec_df)
lmm_varslope_sqrt_log10 = lmer(sqrt(`700`) ~ log10(age) + (log10(age)|scientificName), data = spec_df)
lmm_varslope_log_log10 = lmer(log(`700`) ~ log10(age) + (log10(age)|scientificName), data = spec_df)
lmm_varslope_log10_sqrt = lmer(log10(`700`) ~ sqrt(age) + (sqrt(age)|scientificName), data = spec_df)
lmm_fixedslope = lmer(sqrt(`700`) ~ log10(age) + (1|scientificName), data = spec_df)


lmm_no_log10age = lmer(sqrt(`700`) ~ age + (age|scientificName), data = spec_df)
AIC(lmm_no_log10age)

AIC(linear)
AIC(glm)
AIC(lmm)
AIC(lmm_varslope)
AIC(lmm_varslope_sqrt)
AIC(lmm_varslope_log10)
AIC(lmm_varslope_log)
AIC(lmm_varslope_log10_2)
AIC(lmm_varslope_log_log10)
AIC(lmm_varslope_sqrt_log10) # best out of variable slope lmms
AIC(lmm_fixedslope) #best out of lmms
AIC(lmm_varslope_log10_sqrt)

summary(lmm_varslope_sqrt_log10)
summary(lmm_fixedslope)


(mm_plot <- ggplot(spec_df, aes(x = log10(age), y = sqrt(`700`))) +
        facet_wrap(~scientificName, nrow=5) +   # a panel for each mountain range
        geom_point(alpha = 0.5) +
        theme_classic() +
        geom_line(data = cbind(spec_df, pred = predict(lmm_fixedslope)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
        theme(legend.position = "none",
              panel.spacing = unit(2, "lines"))  # adding space between panels
)

(mm_plot <- ggplot(spec_df, aes(x = log10(age), y = sqrt(`700`))) +
        facet_wrap(~scientificName, nrow=5) +   # a panel for each mountain range
        geom_point(alpha = 0.5) +
        theme_classic() +
        geom_line(data = cbind(spec_df, pred = predict(lmm_varslope_sqrt_log10)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
        theme(legend.position = "none",
              panel.spacing = unit(2, "lines"))  # adding space between panels
)

################################################################################
#Report desired statistics from all wavelengths
################################################################################
library(optimx)
spectra = readRDS('spectra/lichen_spectra.rds')
spectra = spectra[meta(spectra)$age <= 60, ]
#spectra = normalize(spectra)
#4 scans per individual cannot be treated as independent in a linear model, so 
#I'm taking the mean per individual
#spectra = aggregate(spectra, by = meta(spectra)$X, mean, try_keep_txt(mean))
spec_df = as.data.frame(spectra)

varSlope_aic = c()
fixedSlope_aic = c()

for(i in seq(400, 2400, 10)) {
    x = toString(i)
    lmm_varSlope = lmer(spec_df[, x] ~ age + (age|scientificName), data = spec_df, REML = T, 
                        lmerControl(
                            optimizer ='optimx', optCtrl=list(method='nlminb')))
    lmm_fixedSlope = lmer(spec_df[, x] ~ age + (1|scientificName), data = spec_df, REML = T, 
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

ranEffects = as.data.frame(matrix(nrow = 29))
rownames(ranEffects) = sort(unique(spec_df$scientificName))
int_r2_list = c()
slope_r2_list = c()
age_effect_list = c()
age_97.5_list = c()
age_2.5_list = c()


for(i in seq(400, 2400, 1)) {
    x = toString(i)
    lmm = lmer(spec_df[, x] ~ age + (age|scientificName), data = spec_df, REML = T, 
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

saveRDS(stats_list, 'models/lmm_60yrs.rds')

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




################################################################################
#testing
################################################################################

spectra = readRDS('spectra/lichen_spectra.rds')
#spectra = normalize(spectra)
#4 scans per individual cannot be treated as independent in a linear model, so 
#I'm taking the mean per individual
#spectra = aggregate(spectra, by = meta(spectra)$X, mean, try_keep_txt(mean))
spec_df = as.data.frame(spectra)
spec_df$scaled_age = scale(spec_df$age, scale = T, center = T)


spec_scaled = scale(spectra[,400], center = T, scale = T)
lmm1 = lmer(spec_scaled ~ log10(spec_df$age) + (log10(spec_df$age)|spec_df$scientificName))
lmm2 = lmer(sqrt(spec_df[, 400]) ~ log10(spec_df$age) + (log10(spec_df$age)|spec_df$scientificName))
lmm3 = lmer(spec_scaled ~ spec_df$age + (spec_df$age|spec_df$scientificName)) #fails
lmm4 = lmer(spec_scaled ~ spec_df$scaled_age + (spec_df$scaled_age|spec_df$scientificName))
lmm5 = lmer(sqrt(spec_scaled) ~ log10(spec_df$age) + (log10(spec_df$age)|spec_df$scientificName)) #fails
lmm6 = lmer(spec_scaled ~ sqrt(spec_df$age) + (sqrt(spec_df$age)|spec_df$scientificName))
lmm7 = lmer(spec_scaled ~ spec_df$age + (spec_df$age|spec_df$scientificName))


AIC(lmm1)
AIC(lmm2)
AIC(lmm3)
AIC(lmm4)
AIC(lmm5)



ctl = lmerControl(optimizer = 'Nelder_M')


lmm_varslope1 = lmer(`700` ~ age + (age|scientificName), data = spec_df, REML = FALSE, 
                     control = lmerControl(optimizer = "Nelder_Mead"))
lmm_fixedslope = lmer(`700` ~ age + (1|scientificName), data = spec_df, REML = FALSE, 
                      control = lmerControl(optimizer = "Nelder_Mead"))

