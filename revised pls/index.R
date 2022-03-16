################################################################################
# Revised plsda workflow
################################################################################
# Setup
################################################################################

library(parallel)

#functions
runPlsda = readRDS("revised pls/functions/runPlsda.rds")
getNComps = readRDS("revised pls/functions/getNComps.rds")

#spectra
spectra = readRDS("spectra/lichen_spectra.rds")

# set the directories you would like the downsampled and upsampled model 
# objects to be saved to. The upsampled models will be what you draw your
# conclusions from. The downsampled models will be used to determine the optimal
# number of components to use in the upsampling procedure
downSamplingDirectory = "revised pls/data/downsampling"
upSamplingDirectory = "revised pls/data/upsampling"

# choose class you are trying to classify - scientificName, family, order, or class
className = 'scientificName'

# include the effect of age in the model?
includeAge = FALSE

# number of components to start with. I'm arbitrarily choosing 40, but you can
# change this
ncomps = 4 

# choose the number of iterations you want to run plsda for (I suggest 100)
iterations = seq(1, 4)

################################################################################
# Classification - parallelized
################################################################################

# choose number of cores - this code chooses half the cores available on your
# machine
numCores = floor(detectCores()/2)
cluster = makeCluster(numCores)

# this function will run PLSDA and save model objects into the specified directory
# the console will log a list of NULL values, don't worry about this
parLapply(cl = cluster, iterations, runPlsda, spectra = spectra,
          className = className, ncomp = ncomps, resampling = "down",
          include_age = includeAge, directory = downSamplingDirectory)

# get optimal number of components from downsampled model objects
optimalNumComps = getNComps(downSamplingDirectory, ncomps)

# run plsda with upsampling with optimal number of components
parLapply(cl = cluster, iterations, runPlsda, spectra = spectra,
          className = className, ncomp = optimalNumComps, resampling = "up",
          include_age = includeAge, directory = upSamplingDirectory)

