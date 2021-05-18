#add meta data to spectra
################################################################################
#load packages and data
################################################################################
library(spectrolab)
library(rlist)
library(dplyr)

spec_all = readRDS("spectra/spec_all.rds")
data = read.csv("metadata/metadata_joined.csv", stringsAsFactors = F)

################################################################################
#match spectra and meta data names
################################################################################
#add file names to entries with missing file names
broken = data[data$File.Name == "", ]
broken$File.Name = paste(broken$X, broken$year, sep = ".")
data = data[!data$File.Name == "",]
data = rbind(data, broken)

#clean meta data
data = data[!data$Quality == "BAD",]
data = data[!data$Quality == "Not Scanned, bad",]
data$name = paste(data$X, data$year, sep = ".")
rownames(data) = data$name

#clean up spectra names
m = as.matrix(spec_all)
rownames(m) = gsub(".sig", "", rownames(m))
names(spec_all) = rownames(m)
spec_names = as.data.frame(rownames(m))
colnames(spec_names) = "spectra_ID"

#match data and scan names
idx2 = sapply(rownames(data), grep, spec_names$spectra_ID)

#remove data entries that do not have associated spectra
no_spectra = c("654463.1974", "721604.1980", "802685.1988", "815554.1991",
               "22571.1892", "772360.1971", "934442.1993", "691863.1963",
               "941507.1998", "21569.NA", "21435.1897", "673814.1974",
               "15239.1899", "871614.2000", "14293.1899")
data = data[!rownames(data) %in% no_spectra,]

#match names again
idx2 = sapply(rownames(data), grep, spec_names$spectra_ID)
idx1 <- sapply(seq_along(idx2), function(i) rep(i, length(idx2[[i]])))
new_data = cbind(data[unlist(idx1),,drop=F], spec_names[unlist(idx2),,drop=F])

#bring to excell to modify rownames to match spectra names
#write.csv(new_data, "metadata/meta_full.csv")

################################################################################
#Now remove repeated names in scan_names, then remove any names in the metadata 
#that are not in the spectra names
################################################################################
data = read.csv("metadata/meta_full.csv", stringsAsFactors = F)
rownames(data) = data[,1]
data = data[,-c(1,2)]

#remove duplicate spectra names and spectra without metadata
spec_names_no_dup = spec_names %>% distinct()
remove.spec = base::setdiff(spec_names_no_dup$spectra_ID, rownames(data))
spec_names_2 = spec_names_no_dup[!spec_names_no_dup$spectra_ID %in% remove.spec,]

#remove metadata without spectra
remove.data = base::setdiff(rownames(data), spec_names_no_dup$spectra_ID)
data = data[!rownames(data) %in% remove.data, ]

#remove duplicate spectra and spectra without metadata
spec_m = as.matrix(spec_all)
spec_df = as.data.frame(spec_m)
spec_df = spec_df %>% distinct()
rownames(spec_df) = gsub("X", "", rownames(spec_df))
spec_new = as_spectra(spec_df)
remove_spectra = base::setdiff(names(spec_new), spec_names_2)
spec_new = spec_new[!names(spec_new) %in% remove_spectra, ]

#remove metadata that somehow is still different than spectra
remove_meta = base::setdiff(rownames(data), names(spec_new))
data = data[!rownames(data) %in% remove_meta, ]

#set metadata in same order as spectra
data = data[names(spec_new),,drop = F]

#Finally, add the metadata
meta(spec_new) = data

saveRDS(spec_new, "spectra/lichen_spectra.rds")

#fix some typing errors in the metadata and resave
spec_all = readRDS("spectra/lichen_spectra.rds")

meta(spec_all)$Morphology = gsub("Foliose", "foliose", meta(spec_all)$Morphology)
meta(spec_all)$scientificName = gsub(" ","_", meta(spec_all)$scientificName)

saveRDS(spec_all, "spectra/lichen_spectra.rds")
#vector normalize
vn = normalize(spec_all)
saveRDS(vn, "spectra/lichen_spectra_vn.rds")
