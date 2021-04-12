library(dplyr)

################################################################################
#remove bad samples
################################################################################

df = read.csv("metadata/metadata_joined.csv", stringsAsFactors = F)
df = df[!df$File.Name == "",]
df = df[!is.na(df$year),]
df = df[-13,]
df = df[!df$Quality == "BAD",]
df = df[!df$Quality == "not great",]

################################################################################
#create time classes function
################################################################################
#function that creates n time classes and adds a column with those classes
classify_time = function(dataframe, time_name, n_classes, time_span, start_year){
  x = time_span/n_classes
  for (i in 1:n_classes){
    if (i == 1){
      dataframe$timeClass[start_year - dataframe$year <= x] = i
    }
    else if (i > 1) {
      dataframe$timeClass[start_year - dataframe$year > (i-1)*x &
                            start_year - dataframe$year <= i*x] = i
    }
  }
  names(dataframe)[names(dataframe) == "timeClass"] <- time_name
  
  return(dataframe)
}

################################################################################
#create time classes for 126 interval
################################################################################

#time classes, referred to as 'nsplit', from 3 to 12
df = classify_time(dataframe = df, time_name = "split2", n_classes = 2,
                   time_span = 126, start_year = 2012)
df = classify_time(dataframe = df, time_name = "split3", n_classes = 3,
                   time_span = 126, start_year = 2012)
df = classify_time(dataframe = df, time_name = "split4", n_classes = 4,
                   time_span = 126, start_year = 2012)
df = classify_time(dataframe = df, time_name = "split5", n_classes = 5,
                   time_span = 126, start_year = 2012)
df = classify_time(dataframe = df, time_name = "split6", n_classes = 6,
                   time_span = 126, start_year = 2012)
df = classify_time(dataframe = df, time_name = "split7", n_classes = 7,
                   time_span = 126, start_year = 2012)
df = classify_time(dataframe = df, time_name = "split8", n_classes = 8,
                   time_span = 126, start_year = 2012)
df = classify_time(dataframe = df, time_name = "split9", n_classes = 9,
                   time_span = 126, start_year = 2012)
df = classify_time(dataframe = df, time_name = "split10", n_classes = 10,
                   time_span = 126, start_year = 2012)
df = classify_time(dataframe = df, time_name = "split11", n_classes = 11,
                   time_span = 126, start_year = 2012)
df = classify_time(dataframe = df, time_name = "split12", n_classes = 12,
                   time_span = 126, start_year = 2012)


df2 = df %>% group_by(scientificName) %>% filter(length(unique(split2)) == 2)
df3 = df %>% group_by(scientificName) %>% filter(length(unique(split3)) == 3)
df4 = df %>% group_by(scientificName) %>% filter(length(unique(split4)) == 4)
df5 = df %>% group_by(scientificName) %>% filter(length(unique(split5)) == 5)
df6 = df %>% group_by(scientificName) %>% filter(length(unique(split6)) == 6)
df7 = df %>% group_by(scientificName) %>% filter(length(unique(split7)) == 7)
df8 = df %>% group_by(scientificName) %>% filter(length(unique(split8)) == 8)

species_10 = c('Caloplaca flavovirescens', 'Candelaria concolor',
               'Dimelaena oreina', 'Flavoparmelia caperata', 'Graphis scripta',
               'Ionaspis lacustris', 'Peltigera elisabethae',
               'Trypethelium virens', 'Umbilicaria muehlenbergii', 
               'Verrucaria fuscella')

#best two datasets to evaluate the 126 year interval

#course measurement of change with time but with broad taxonomic coverage
df10 = df[df$scientificName %in% species_10,]
write.csv(df10, "metadata/ten_species_metadata.csv")

#Fine measurement of change over 126 years with only one species
tv = df[df$scientificName == "Trypethelium virens",]
write.csv(tv, "metadata/T_virens_metadata.csv")

################################################################################
#time classes 1970-2012
################################################################################
df.new = df[df$year > 1970,]

df.new = classify_time(df.new, time_name = "split10", n_classes = 10,
                         time_span = 42, start_year = 2012)
df.new = classify_time(df.new, time_name = "split9", n_classes = 9,
                       time_span = 42, start_year = 2012)
df.new = classify_time(df.new, time_name = "split8", n_classes = 8,
                       time_span = 42, start_year = 2012)
df.new = classify_time(df.new, time_name = "split7", n_classes = 7,
                       time_span = 42, start_year = 2012)
df.new = classify_time(df.new, time_name = "split6", n_classes = 6,
                       time_span = 42, start_year = 2012)
df.new = classify_time(df.new, time_name = "split5", n_classes = 5,
                       time_span = 42, start_year = 2012)
df.new = classify_time(df.new, time_name = "split4", n_classes = 4,
                       time_span = 42, start_year = 2012)
df.new = classify_time(df.new, time_name = "split3", n_classes = 3,
                       time_span = 42, start_year = 2012)


df.new10 = df.new %>% group_by(scientificName) %>%
  filter(length(unique(split10)) == 10)

df.new9 = df.new %>% group_by(scientificName) %>%
  filter(length(unique(split9)) == 9)

df.new8 = df.new %>% group_by(scientificName) %>%
  filter(length(unique(split8)) == 8)

df.new7 = df.new %>% group_by(scientificName) %>%
  filter(length(unique(split7)) == 7)

df.new6 = df.new %>% group_by(scientificName) %>%
  filter(length(unique(split6)) == 6)

df.new5 = df.new %>% group_by(scientificName) %>%
  filter(length(unique(split5)) == 5)

df.new4 = df.new %>% group_by(scientificName) %>%
  filter(length(unique(split4)) == 4)

df.new3 = df.new %>% group_by(scientificName) %>%
  filter(length(unique(split3)) == 3)

split_3_1 = as.data.frame(table(df.new3[df.new3$split3 == 1,]$scientificName))
split_3_1 = split_3_1[split_3_1$Freq > 2,]
a = as.character(split_3_1$Var1)

split_3_2 = as.data.frame(table(df.new3[df.new3$split3 == 2,]$scientificName))
split_3_2 = split_3_2[split_3_2$Freq > 2,]
b = as.character(split_3_2$Var1)

split_3_3 = as.data.frame(table(df.new3[df.new3$split3 == 3,]$scientificName))
split_3_3 = split_3_3[split_3_3$Freq > 2,]
c = as.character(split_3_3$Var1)

a = intersect(a,b)
b = intersect(b,a)
c = intersect(c,a)

df_42 = df[df$scientificName %in% c,]
write.csv(df_42, "metadata/42_years_3_classes.csv")

################################################################################
#time classes 1990-2010
################################################################################
df.young = df[df$year > 1990,]

df.young = classify_time(df.young, time_name = "split10", n_classes = 10,
                       time_span = 20, start_year = 2010)
df.young10 = df.young %>% group_by(scientificName) %>%
  filter(length(unique(split10)) == 10)

df.young = classify_time(df.young, time_name = "split9", n_classes = 9,
                         time_span = 20, start_year = 2010)
df.young9 = df.young %>% group_by(scientificName) %>%
  filter(length(unique(split9)) == 9)

df.young = classify_time(df.young, time_name = "split8", n_classes = 8,
                         time_span = 20, start_year = 2010)
df.young8 = df.young %>% group_by(scientificName) %>%
  filter(length(unique(split8)) == 8)

df.young = classify_time(df.young, time_name = "split7", n_classes = 7,
                         time_span = 20, start_year = 2010)
df.young7 = df.young %>% group_by(scientificName) %>%
  filter(length(unique(split7)) == 7)

df.young = classify_time(df.young, time_name = "split6", n_classes = 6,
                         time_span = 20, start_year = 2010)
df.young6 = df.young %>% group_by(scientificName) %>%
  filter(length(unique(split6)) == 6)

df.young = classify_time(df.young, time_name = "split5", n_classes = 5,
                         time_span = 20, start_year = 2010)
df.young5 = df.young %>% group_by(scientificName) %>%
  filter(length(unique(split5)) == 5)

df.young = classify_time(df.young, time_name = "split4", n_classes = 4,
                         time_span = 20, start_year = 2010)
df.young4 = df.young %>% group_by(scientificName) %>%
  filter(length(unique(split4)) == 4)

df.young = classify_time(df.young, time_name = "split3", n_classes = 3,
                         time_span = 20, start_year = 2010)
df.young3 = df.young %>% group_by(scientificName) %>%
  filter(length(unique(split3)) == 3)

df.young = classify_time(df.young, time_name = "split2", n_classes = 2,
                         time_span = 20, start_year = 2010)
df.young2 = df.young %>% group_by(scientificName) %>%
  filter(length(unique(split2)) == 2)

split_2_1 = as.data.frame(table(df.young2[df.young2$split2 == 1,]$scientificName))
split_2_1 = split_2_1[split_2_1$Freq > 2,]
a = as.character(split_2_1$Var1)

split_2_2 = as.data.frame(table(df.young2[df.young2$split2 == 2,]$scientificName))
split_2_2 = split_2_2[split_2_2$Freq > 2,]
b = as.character(split_2_2$Var1)

a = intersect(a,b)

df_20 = df[df$scientificName %in% a,]
write.csv(df_20, "metadata/20_years_2_classes.csv")

################################################################################
#time classes 2000-2010
################################################################################
df.x = df[df$year > 2000,]

df.x = classify_time(df.x, time_name = "split10", n_classes = 10,
                         time_span = 10, start_year = 2010)
df.x10 = df.x %>% group_by(scientificName) %>%
  filter(length(unique(split10)) == 10)

df.x = classify_time(df.x, time_name = "split9", n_classes = 9,
                         time_span = 10, start_year = 2010)
df.x9 = df.x %>% group_by(scientificName) %>%
  filter(length(unique(split9)) == 9)

df.x = classify_time(df.x, time_name = "split8", n_classes = 8,
                         time_span = 10, start_year = 2010)
df.x8 = df.x %>% group_by(scientificName) %>%
  filter(length(unique(split8)) == 8)

df.x = classify_time(df.x, time_name = "split7", n_classes = 7,
                         time_span = 10, start_year = 2010)
df.x7 = df.x %>% group_by(scientificName) %>%
  filter(length(unique(split7)) == 7)

df.x = classify_time(df.x, time_name = "split6", n_classes = 6,
                         time_span = 10, start_year = 2010)
df.x6 = df.x %>% group_by(scientificName) %>%
  filter(length(unique(split6)) == 6)

df.x = classify_time(df.x, time_name = "split5", n_classes = 5,
                         time_span = 10, start_year = 2010)
df.x5 = df.x %>% group_by(scientificName) %>%
  filter(length(unique(split5)) == 5)

df.x = classify_time(df.x, time_name = "split4", n_classes = 4,
                     time_span = 10, start_year = 2010)
df.x4 = df.x %>% group_by(scientificName) %>%
  filter(length(unique(split4)) == 4)

df.x = classify_time(df.x, time_name = "split3", n_classes = 3,
                     time_span = 10, start_year = 2010)
df.x3 = df.x %>% group_by(scientificName) %>%
  filter(length(unique(split3)) == 3)

df.x = classify_time(df.x, time_name = "split2", n_classes = 2,
                     time_span = 10, start_year = 2010)
df.x2 = df.x %>% group_by(scientificName) %>%
  filter(length(unique(split2)) == 2)
