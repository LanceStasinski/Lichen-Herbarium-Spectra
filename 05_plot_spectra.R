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
par(mfrow = c(1,1))
spec15 = spec[meta(spec)$age <= 15, ]
spec30_1 = spec[meta(spec)$age > 15, ] 
spec30 = spec30_1[meta(spec30_1)$age <= 30,]
spec45_1 = spec[meta(spec)$age > 30, ]
spec45 = spec45_1[meta(spec45_1)$age <= 45,]
spec60_1 = spec[meta(spec)$age > 45, ]
spec60 = spec60_1[meta(spec60_1)$age <= 60,]

jpeg(filename = '../../lichen figures/age_spec_mean.jpeg',
     width = 6, height = 4, units = 'in', res = 1200)
plot(mean(spec15), col = '#018571', ylab = 'Reflectance', xlab = 'Wavelength (nm)')
plot(mean(spec30), col = '#80cdc1', add = T)
plot(mean(spec45), col = '#dfc27d', add = T)
plot(mean(spec60), col = '#a6611a', add = T)
dev.off()

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
par(mfrow = c(2,2))
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

#Caloplaca flavovirescens
spec_c = spec[meta(spec)$scientificName == 'Caloplaca_flavovirescens',]
spec_c = spec_c[meta(spec_c)$age < 60, ]

spec15 = spec_c[meta(spec_c)$age <= 15, ]
spec30_1 = spec_c[meta(spec_c)$age > 15, ] 
spec30 = spec30_1[meta(spec30_1)$age <= 30,]
spec45_1 = spec_c[meta(spec_c)$age > 30, ]
spec45 = spec45_1[meta(spec45_1)$age <= 45,]
spec60_1 = spec_c[meta(spec_c)$age > 45, ]
spec60 = spec60_1[meta(spec60_1)$age <= 60,]


plot(mean(spec15), col = 'green')
plot(mean(spec30), col = 'darkgreen', add = T)
plot(mean(spec45), col = 'brown', add = T)
plot(mean(spec60), col = 'gray', add = T)

#Loxospora_elatina
spec_l = spec[meta(spec)$scientificName == 'Loxospora_elatina',]
spec_l = spec_l[meta(spec_l)$age < 60, ]

spec15 = spec_l[meta(spec_l)$age <= 15, ]
spec30_1 = spec_l[meta(spec_l)$age > 15, ] 
spec30 = spec30_1[meta(spec30_1)$age <= 30,]
spec45_1 = spec_l[meta(spec_l)$age > 30, ]
spec45 = spec45_1[meta(spec45_1)$age <= 45,]
spec60_1 = spec_l[meta(spec_l)$age > 45, ]
spec60 = spec60_1[meta(spec60_1)$age <= 60,]
spec60 = NULL

par(mfrow = c(1,1))
plot(mean(spec15), col = 'green')
plot(mean(spec30), col = 'darkgreen', add = T)
plot(mean(spec45), col = 'brown', add = T)
if (!is.null(spec60)){
  plot(mean(spec60), col = 'gray', add = T)
}


#Chrysothrix_candelaris
spec_ch = spec[meta(spec)$scientificName == 'Chrysothrix_candelaris',]
spec_ch = spec_ch[meta(spec_ch)$age < 60, ]

spec15 = spec_ch[meta(spec_ch)$age <= 15, ]
spec30_1 = spec_ch[meta(spec_ch)$age > 15, ] 
spec30 = spec30_1[meta(spec30_1)$age <= 30,]
spec45_1 = spec_ch[meta(spec_ch)$age > 30, ]
spec45 = spec45_1[meta(spec45_1)$age <= 45,]
spec60_1 = spec_ch[meta(spec_ch)$age > 45, ]
spec60 = spec60_1[meta(spec60_1)$age <= 60,]


plot(mean(spec15), col = 'green')
plot(mean(spec30), col = 'darkgreen', add = T)
plot(mean(spec45), col = 'brown', add = T)
plot(mean(spec60), col = 'gray', add = T)

#all species over time
#function
plot_species_spec = function(species_spec, age) {
  if (age == 15) {
    spec15 = NULL
    spec15 = species_spec[meta(species_spec)$age <= 15, ]
    spec
    if (!is.null(spec15)) {
      spec15 = aggregate(spec15, by = meta(spec15)$X, mean, try_keep_txt(mean))
      plot(spec15, col = '#018571', add = T)
    }
  } else if (age == 30) {
    spec30_1 = NULL
    spec30 = NULL
    spec30_1 = species_spec[meta(species_spec)$age > 15, ]
    spec30 = spec30_1[meta(spec30_1)$age <= 30,]
    
    if (!is.null(spec30)) {
      spec30 = aggregate(spec30, by = meta(spec30)$X, mean, try_keep_txt(mean))
      plot(spec30, col = '#80cdc1', add = T)
    }
  } else if (age == 45) {
    spec45_1 = NULL
    spec45 = NULL
    spec45_1 = species_spec[meta(species_spec)$age > 30, ]
    spec45 = spec45_1[meta(spec45_1)$age <= 45,]
    if (!is.null(spec45)) {
      spec45 = aggregate(spec45, by = meta(spec45)$X, mean, try_keep_txt(mean))
      plot(spec45, col = '#dfc27d', add = T)
    }
  } else if (age == 60) {
    spec60_1 = NULL
    spec60 = NULL
    spec60_1 = species_spec[meta(species_spec)$age > 45, ]
    spec60 = spec60_1[meta(spec60_1)$age <= 60,]
    if (!is.null(spec60)) {
      spec60 = aggregate(spec60, by = meta(spec60)$X, mean, try_keep_txt(mean))
      plot(spec60, col = '#a6611a', add = T)
    }
  }
  
}


plot_aging_spec = function(species) {
  
  for(i in 1:length(species)) {
    species_spec = spec[meta(spec)$scientificName == species[i],]

    null_spec = rep(-.05, 2001)
    m = matrix(nrow = 1, ncol = 2001)
    colnames(m) = seq(400, 2400, 1)
    m = rbind(m, null_spec)
    m = as.matrix(m[-1,])
    null_spec = as_spectra(t(m))
    plot(null_spec, ylim = c(0, 0.8), xlab = 'Wavelength (nm)',
         ylab = 'Reflectance', main = species[i])
    
    try(plot_species_spec(species_spec, 15))
    try(plot_species_spec(species_spec, 30))
    try(plot_species_spec(species_spec, 45))
    try(plot_species_spec(species_spec, 60))
    
  }
}




#plots
spec = readRDS("spectra/lichen_spectra.rds")
spec = spec[meta(spec)$age <= 60,]
species = sort(unique(meta(spec)$scientificName))

par(mfrow = c(3,2))

plot_aging_spec(species)







