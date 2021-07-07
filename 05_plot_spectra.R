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

