## define the first lag autocorrelation coefficient to prewhite the data as an AR(Kmax)function
## 
##  Input:
##     data: one dimensional daily timeseries
##   
##  output:
##    ak.lag: fist lag autocorrelation coefficient
##    data.prewhite: matrix after removing of the first lag autorcorrelation
##                   if ak_lag is ss, original data otherwise.
##    ak.ss (integer): statistical significance of the first lag autocorrelation:
##                     alpha.ak is ss at the alpha.ak level, zero otherwise

## Martine Collaud Coen, MeteoSwiss, 9.2019
## Alessandro Bigi, University of Modena and Reggio Emilia, Sep 2020
## data input test is  data  <-  unlist(read.table("../test_data/nanprewhite_AR_test1_in.csv"))
## test on Oct 12th: ok!


nanprewhite.AR <- function(data, alpha.ak = 95){

    require(signal)

    if (alpha.ak > 100) stop('the confidence limit has to be lower than 100%')

    
    if (sum(is.na(data)) == length(data)){

        ak.lag <- NA
        data.prewhite <- NA
        
    } else {

        ## remove infinite data if existing
        data[which(is.infinite(data) == TRUE)] <- NA
        data[which(is.nan(data) == TRUE)] <- NA
        
        ## specifies the number of lags supposed to have a significant autocorrelation coefficient
        p <- 5
        
        ##  number of lags to be computed
        nblag <- 10
        
        ## restrict the number of computed lag and the number of significant AC coef. if the time series is short
        if (nblag > length(data) / 2){
            nblag  <- floor( length(data) / 2 )
        }
        
        if (p >= nblag) {
            p <- nblag-1
        }

        x <- nanautocorr(data, nblag, p)
                
        x1  <-  x / length(na.omit(data))
        lev  <-  signal::levinson(x1, p)
        K <- lev$ref
        
        ak.coef  <- -K
        uconf <- qnorm( 1 - (1-alpha.ak / 100) /2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) / sqrt(sum(complete.cases(data)))

        ak.lag <- x[2]

        ## compute the prewhitened data only is ak(1) is s.s. at 95% confidence limit 
        ## other wise the orgininal is given as an output
        if (abs( ak.coef[1] ) - uconf < 0){
            data.prewhite <- data
            ak.ss <- 0 
            warning('no s.s. autocorrelation')
        } else {
            ##return(list(data, x))
            y <- rep(NA, length(data))
            y[-1] <- x[2] * data[ -length(data) ]
            data.prewhite <- data - y
            ak.ss <- alpha.ak
        }
        
        return(list("ak.lag" = ak.lag, "ak.ss" = ak.ss, "data.prewhite" = data.prewhite))
    }
}
