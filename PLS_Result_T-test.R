################################################################################
#Determine if age makes a significant difference in model accuracy - t-test
################################################################################

# function: from https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha
# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}

#################
#species
#################

t.test2(m1 = 76.94587, m2 = 80.49003, s1 = 1.97579, s2 = 1.91673, n1 = 100, n2 = 100)

#################
#family 
#################

t.test2(m1 = 94.50704, m2 = 94.43662, s1 = 1.052847, s2 = 1.252978, n1 = 100, n2 = 100)

#################
#order
#################

t.test2(m1 = 89.0056, m2 = 89.493, s1 = 1.28294, s2 = 1.380609, n1 = 100, n2 = 100)

#################
#class
#################

t.test2(m1 = 82.45833, m2 = 83.86111, s1 = 2.141263, s2 = 2.012594, n1 = 100, n2 = 100)





