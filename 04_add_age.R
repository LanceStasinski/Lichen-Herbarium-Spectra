#Add age column to metadata
################################################################################
#setup
################################################################################

library(spectrolab)
spectra = readRDS('spectra/lichen_spectra.rds')

################################################################################
#add age
################################################################################
data = meta(spectra)
data$age = 2019 - data$year

m = as.matrix(spectra)
spectra.new = as_spectra(m)
meta(spectra.new) = data

saveRDS(spectra.new, 'spectra/lichen_spectra.rds')