################################################################################
#read in data
################################################################################

taxonomy = read.csv("metadata/taxonomy.csv", stringsAsFactors = F)
spectra_metadata = read.csv("metadata/spectra_metadata.csv",
                            stringsAsFactors = F)

################################################################################
#make sure species names match
################################################################################
taxonomy$Species = gsub("Chysothrix candelaris", "Chrysothrix candelaris",
                        taxonomy$Species)
taxonomy$Species = gsub("Flavopunctelia darrowi",
                                       "Flavopunctelia darrowii", 
                                       taxonomy$Species)

notInMeta = setdiff(unique(taxonomy$Species),
                    unique(spectra_metadata$scientificName))
notInTaxa = setdiff(unique(spectra_metadata$scientificName),
                    unique(taxonomy$Species))

names(taxonomy)[names(taxonomy) == "Species"] <- "scientificName"
################################################################################
#join data
################################################################################

metadata = merge(spectra_metadata, taxonomy, by = "scientificName")

#remove incomplete taxonomy columns
metadata = subset(metadata, select = -c(class, order, family, Chemistry.Methods,
                                        endDayOfYear))
write.csv(metadata, "metadata/metadata_joined.csv")
