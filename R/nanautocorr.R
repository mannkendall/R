#' Autocorrelation function
#'
#' `nanautocorr` computes autocorrelation coefficient dealing with missing values
#' 
#' Calculates the `nlag` autocorrelation coefficient for a data vector containing missing values
#' couples of data including `NA` are excluded from the computation.
#' Here the ACF is calculated using the Pearson\'s correlation coefficient for each lag.
#' 
#' @param data vector of time series
#' @param nlags maximum number of lags in the computation of ACF
#' @param R number of lags until the model is  supposed to have a significant AC coefficient (optional)
#' @return out vector of the ACF
#' @author Fabio Oriani, fabio.oriani@unine.ch, University of Neuchatel (CH) original author; Alessandro Bigi (abigi@unimore.it), University of Modena and Reggio Emilia (IT) for the R version
#'
#' @export

nanautocorr <- function(data, nlags, R){

    if (sum(is.na((data))) > length(data) / 3) {
        
        warning("AC: too many nans","more than a third of data is NA. Autocorrelation is not reliable")
    }
    
    out <- rep(NA, nlags)
    out[1] <- 1
    data <- data - mean(data, na.rm = TRUE)

    ## use segnan to make several input for lpc
    ## without NaNs
   
    for (i in 2 : c(nlags+1)){
        out[i] <- cor(data[i:length(data)], data[1:c(length(data) - i + 1)], use="complete")
    }

    nargin <- length(as.list(match.call())) - 1 ## number of arguments provided
    
    if (nargin == 3){
        if (R >= nlags){
            error('R must be lower than nlags')
        }
    }
    if (missing(R)) R <- 0
        
    ## confidence bounds
    b <- 1.96 * length(data)^(-0.5) * sum(out[1:c(R + 1)]^2) ^(0.5)

    out
}
