################################################################################
# Revised plsda workflow
################################################################################
# Setup
################################################################################

library(parallel)

runPlsda = readRDS("revised pls/functions/runPlsda.rds")
spectra = readRDS('spectra/lichen_spectra.rds')

################################################################################
# Classification - parallelized
################################################################################

# choose number of cores - this code chooses half the cores available on your
# machine
numCores = detectCores()
cluster = makeCluster(names = floor(numCores/2))

# choose the number of iterations you want (I suggest 100)
iterations = seq(1, 4)

# this function will run PLSDA and save model objects into the specified directory
# the console will log a list of NULL values, don't worry about this
parLapply(cl = cluster, iterations, runPlsda, spectra = spectra,
          className = "scientificName", ncomp = 3, resampling = 'down',
          include_age = F, directory = "revised pls/data")
