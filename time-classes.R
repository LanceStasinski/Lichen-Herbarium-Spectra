################################################################################
#remove bad samples
################################################################################

df = read.csv("metadata/spectra_metadata.csv", stringsAsFactors = F)
df = df[!df$File.Name == "",]
df = df[!is.na(df$year),]
df = df[-13,]
df = df[!df$Quality == "BAD",]
df = df[!df$Quality == "not great",]

################################################################################
#create time classes
################################################################################
#6 time classes over 126 years
df["period1"] <- NA
x = 126/6
df$period1[2012 - df$year <= x] = 1
df$period1[2012 - df$year > x & 2012 - df$year <= 2*x] = 2
df$period1[2012 - df$year > 2*x & 2012 - df$year <= 3*x] = 3
df$period1[2012 - df$year > 3*x & 2012 - df$year <= 4*x] = 4
df$period1[2012 - df$year > 4*x & 2012 - df$year <= 5*x] = 5
df$period1[2012 - df$year > 5*x] = 6

table(df$period1)

#5 time classes over 126 years
df["period2"] <- NA

df$period2[2012 - df$year <= x] = 1
df$period2[2012 - df$year > x & 2012 - df$year <= 2*x] = 2
df$period2[2012 - df$year > 2*x & 2012 - df$year <= 3*x] = 3
df$period2[2012 - df$year > 3*x & 2012 - df$year <= 4*x] = 4
df$period2[2012 - df$year > 4*x] = 5
table(df$period2)

#4 time classes over 126 years

df["period3"] <- NA

df$period3[2012 - df$year <= 10] = 1
df$period3[2012 - df$year > x & 2012 - df$year <= 2*x] = 2
df$period3[2012 - df$year > 2*x & 2012 - df$year <= 3*x] = 3
df$period3[2012 - df$year > 3*x] = 4
table(df$period3)


classify_time = function(dataframe, time_name, n_classes, time_span, start_year){
  dataframe[time_name] = NA
  x = time_span/n_classes
  for (i in 1:n_classes){
    if (i == 1){
      dataframe$time_name[start_year - dataframe$year <= x] = i
    }
    else if (i > 1) {
      dataframe$time_name[start_year - dataframe$year > (i-1)*x &
                            start_year - dataframe$year <= i*x] = i
    }
  }
}

classify_time(dataframe = df, time_name = "period1", n_classes = 5, time_span = 126, start_year = 2012)






