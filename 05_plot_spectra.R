#plot spectra
################################################################################
#load packages and data
################################################################################
library(spectrolab)

spec = readRDS("spectra/lichen_spectra.rds")
vn_spec = readRDS("spectra/lichen_spectra_vn.rds")

################################################################################
#spec by species
################################################################################

plot_spectra = function(spectra) {
  species = unique(meta(spectra)$scientificName)
  
  par(mfrow = c(5,6))
  for (i in 1:length(species)){
    plot(spectra[meta(spectra)$scientificName == species[i],],
         xlab = "Wavelength (nm)", 
         ylab = "Reflectance",
         main = species[i])
  }
}

plot_spectra(spectra = spec)
plot_spectra(spectra = vn_spec)

################################################################################
#spec by growth form
################################################################################
meta(spec)$Morphology = gsub("Foliose", "foliose", meta(spec)$Morphology)
meta(vn_spec)$Morphology = gsub("Foliose", "foliose", meta(vn_spec)$Morphology)


plot_spectra = function(spectra) {
  growth = unique(meta(spectra)$Morphology)
  
  par(mfrow = c(2,2))
  for (i in 1:length(growth)){
    plot(spectra[meta(spectra)$Morphology == growth[i],],
         xlab = "Wavelength (nm)", 
         ylab = "Reflectance",
         main = growth[i])
  }
}

plot_spectra(spectra = spec)
plot_spectra(spectra = vn_spec)

################################################################################
#spec by order
################################################################################

plot_spectra = function(spectra) {
  order = unique(meta(spectra)$Order)
  
  par(mfrow = c(4,4))
  for (i in 1:length(order)){
    plot(spectra[meta(spectra)$Order == order[i],],
         xlab = "Wavelength (nm)", 
         ylab = "Reflectance",
         main = order[i])
  }
}

plot_spectra(spectra = spec)
plot_spectra(spectra = vn_spec)

################################################################################
#spec by class
################################################################################

plot_spectra = function(spectra) {
  class = unique(meta(spectra)$Class)
  
  par(mfrow = c(3,2))
  for (i in 1:length(class)){
    plot(spectra[meta(spectra)$Class == class[i],],
         xlab = "Wavelength (nm)", 
         ylab = "Reflectance",
         main = class[i])
  }
}

plot_spectra(spectra = spec)
plot_spectra(spectra = vn_spec)

################################################################################
#spec by age
################################################################################
# 10 years
spec10 = spec[meta(spec)$age <= 10, ]
spec20_1 = spec[meta(spec)$age > 10, ] 
spec20 = spec20_1[meta(spec20_1)$age <= 20,]
spec30_1 = spec[meta(spec)$age > 20, ]
spec30 = spec30_1[meta(spec30_1)$age <= 30,]
spec40_1 = spec[meta(spec)$age > 30, ]
spec40 = spec40_1[meta(spec40_1)$age <= 40,]
spec50_1 = spec[meta(spec)$age > 40, ]
spec50 = spec50_1[meta(spec50_1)$age <= 50,]
spec60_1 = spec[meta(spec)$age > 50, ]
spec60 = spec60_1[meta(spec60_1)$age <= 60,]

plot(mean(spec10), col = 'green')
plot(mean(spec20), col = 'darkgreen', add = T)
plot(mean(spec30), col = 'darkseagreen', add = T)
plot(mean(spec40), col = 'burlywood', add = T)
plot(mean(spec50), col = 'brown', add = T)
plot(mean(spec60), col = 'gray', add = T)

# 15 years
spec15 = spec[meta(spec)$age <= 15, ]
spec30_1 = spec[meta(spec)$age > 15, ] 
spec30 = spec30_1[meta(spec30_1)$age <= 30,]
spec45_1 = spec[meta(spec)$age > 30, ]
spec45 = spec45_1[meta(spec45_1)$age <= 45,]
spec60_1 = spec[meta(spec)$age > 45, ]
spec60 = spec60_1[meta(spec60_1)$age <= 60,]


plot(mean(spec15), col = 'green')
plot(mean(spec30), col = 'darkgreen', add = T)
plot(mean(spec45), col = 'brown', add = T)
plot(mean(spec60), col = 'gray', add = T)

# 20 years
spec20 = spec[meta(spec)$age <= 20, ]
spec40_1 = spec[meta(spec)$age > 20, ] 
spec40 = spec40_1[meta(spec40_1)$age <= 40,]
spec60_1 = spec[meta(spec)$age > 20, ]
spec60 = spec60_1[meta(spec60_1)$age <= 60,]


plot(mean(spec20), col = 'green')
plot(mean(spec40), col = 'brown', add = T)
plot(mean(spec60), col = 'gray', add = T)

# Flavoparmelia caperata
spec_f = spec[meta(spec)$scientificName == 'Flavoparmelia_caperata',]
spec_f = spec_f[meta(spec_f)$age < 60, ]

spec15 = spec_f[meta(spec_f)$age <= 15, ]
spec30_1 = spec_f[meta(spec_f)$age > 15, ] 
spec30 = spec30_1[meta(spec30_1)$age <= 30,]
spec45_1 = spec_f[meta(spec_f)$age > 30, ]
spec45 = spec45_1[meta(spec45_1)$age <= 45,]
spec60_1 = spec_f[meta(spec_f)$age > 45, ]
spec60 = spec60_1[meta(spec60_1)$age <= 60,]


plot(mean(spec15), col = 'green')
plot(mean(spec30), col = 'darkgreen', add = T)
plot(mean(spec45), col = 'brown', add = T)
plot(mean(spec60), col = 'gray', add = T)
