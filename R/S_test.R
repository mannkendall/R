## compute the S statistic (Si) for the MK test

## INPUT
## "data" is the vector of observations
## "t.time" is a matrix with the year, month, date vector as individual variables in each column. The number of rows is the number of obsverations

##  "nb.an", "an.debut", "an.end" are respectively the number of years, the start year and the end year of the series

## OUTPUT
## "S" is the double sum on the sign of the difference between data (Si)
## "n" is the number of valid data in each year of the time series

## R is the "rank". Without ties and NaN, R is the Spearman correlation coefficient

## source: Gilbert, 1987

## OK test Oct 16

S.test <- function(data, t.time){

    ## if length(data)~=length(time)
    ##     error(' data and time should have the same length');
    ## end
    
    ## determine the non-missing data number
    L <- length(data)

    ## start year
    an.debut <- min(t.time[,1], na.rm = TRUE)
    ## end year
    an.end  <- max(t.time[,1], na.rm = TRUE)
    ## number of years
    nb.an <- an.end - an.debut + 1
 
    ## initialise matrix Sij
    Sij <- matrix(NA, nrow = nb.an, ncol = nb.an)
    n <- rep(NA, nb.an)

    for (i in an.debut:an.end) {

        ## full length of the valid data in the series for the i-th year
        n[i - an.debut + 1] <- length( na.omit(data[t.time[,1]==i]) )

        if (i < an.end){

            for (j in c(i + 1) : an.end) {

                nj <- length(data[t.time[,1]==j])
                sj <- rep(NA, nj)
            
                dataj <- data[t.time[,1]==j]

                for (k in 1:nj){
                    sj[k] <- sum( sign( dataj[k] - data[t.time[,1]==i]), na.rm=TRUE)
                }
                
                Sij[j - an.debut + 1, i - an.debut + 1] <- sum(sj, na.rm=TRUE)
            }
        }
    }
    
    S <- sum(Sij, na.rm = TRUE)   
    return(list("n"=n, "S"=S))
}
