#Process the spectra: add metadata, trim spectra
################################################################################
#load packages
################################################################################

library(spectrolab)
library(rlist)

################################################################################
#clean spectra fucntion
################################################################################
clean_spectra = function(spectra) {
  #match sensor overlap
  matched = match_sensors(spectra, splice_at = c(990, 1900),
                          interpolate_wvl = 10)
  #trim spectra to 400:2400 nm
  trimmed = matched[ , bands(matched, 400, 2400)]
  
  #resample to give all spectra same number of bands. 2.231 chosen for band 
  #resolution because this is the average band size of 896 bands between 
  #400.7 and 2399.4. 896 is the lowest band count for these spectra, and the 
  #different spectra within this dataset may not be the same size due to a 
  #change in collection procedure or change in the selected parameters of on the
  #instrument
  resampled = resample(trimmed, seq(400.7, 2399.4, 2.231))
  
  #remove any unlabeled white references
  noWR = resampled[!rowSums(resampled > 1),]
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
spec_all = Reduce(combine, spec_list)

#handle the peltigera spectra
spec24 = read_spectra(path = spec.dirs[24], format = "sig", 
                      exclude_if_matches = c("BAD", "WR"))
spec24.1 = spec24[[1]]
spec24.2 = spec24[[2]]

spec24.1.1 = clean_spectra(spec24.1)
spec24.2.1 = clean_spectra(spec24.2)

#add peltigera spectra to full spectra
spec_all = Reduce(combine, list(spec_all, spec24.1.1, spec24.2.1))

#smooth
spec_all = smooth(spec_all)

vn = normalize(spec_all)








