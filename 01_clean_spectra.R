#Process the spectra: add metadata, trim spectra
################################################################################
#load packages
################################################################################

library(spectrolab)
library(rlist)

################################################################################
#clean spectra function
################################################################################
clean_spectra = function(spectra) {
  #match sensor overlap if overlapped regions occur
  
  if (ncol(as.matrix(spectra)) == 1024) {
  matched = match_sensors(spectra, splice_at = c(990, 1900),
                          interpolate_wvl = c(5,1))
  
  #Resample to 1nm resolution to make interpretation easier and trim to 400:2400nm
  trimmed = matched[, bands(matched, 392, 2415)]
  resampled = resample(trimmed, seq(392,2415, 1))
  t.resamp = resampled[, 400:2400]
  } else {
    trimmed = spectra[, bands(spectra, 392, 2415)]
    resampled = resample(trimmed, seq(392,2415, 1))
    t.resamp = resampled[, 400:2400]
  }
  
  #remove any unlabeled white references
  noWR = t.resamp[!rowSums(t.resamp > 1),]
  return(noWR)
}

################################################################################
#read in and clean spectra
################################################################################
spec.dirs = list.dirs(path = "./spectra")
spec.dirs = spec.dirs[-1]
spec.dirs2 = spec.dirs[-24] #peltigera elisabethae has two sets of spectra:
#1 with 1024 bands and another with 996 bands. These spectra will be handled 
#separately from the rest

spec_list = list()
for (i in 1:length(spec.dirs2)) {
  #read in spectra. Exclude bad and white reference scans
  raw = read_spectra(path = spec.dirs2[i], format = "sig", 
                      exclude_if_matches = c("BAD", "WR"))
  
  #match sensors, trim and resample spectra
  clean_spec = clean_spectra(raw)
  
  #add spectra to list of spectra
  spec = assign(paste0("spec", i), clean_spec)
  spec_list = list.append(spec_list, get("spec"))
}
#combine spectra into single spectra object
spec_all = Reduce(spectrolab::combine, spec_list)

#handle the peltigera spectra
spec.p = read_spectra(path = spec.dirs[24], format = "sig", 
                      exclude_if_matches = c("BAD", "WR"))
spec.p.1 = spec.p[[1]]
spec.p.2 = spec.p[[2]]

spec.p.1.1 = clean_spectra(spec.p.1)

spec.p.2.1 = spec.p.2[, bands(spec.p.2, 392, 2415)]
spec.p.2.2 = resample(spec.p.2.1, seq(392,2415, 1))
spec.p.2.3 = spec.p.2.2[, 400:2400]


#add peltigera spectra to full spectra
spec_all = Reduce(spectrolab::combine, list(spec_all, spec.p.1.1, spec.p.2.3))

#smooth
spec_all = smooth(spec_all)
saveRDS(spec_all, "spectra/spec_all.rds")










