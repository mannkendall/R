#' Compute AR1 prior to prewhite
#'
#' Define the first lag autocorrelation coefficient to prewhite the data as an AR(Kmax) function
#' 
#' @param data vector to be analysed
#' @param alpha.ak statistical significance in percentage for the first lag autocorrelation (default is 95)  
#'
#' @return ak.lag fist lag autocorrelation coefficient
#' @return data.prewhite data after removing of the first lag autorcorrelation; if ak.lag is statistically significant, original data otherwise
#' @return ak.ss statistical significance of the first lag autocorrelation: alpha.ak is the statistical significance at the alpha.ak level, zero otherwise
#'
#' @author Martine Collaud Coen (martine.collaud@meteoswiss.ch), MeteoSwiss (CH) and Alessandro Bigi (abigi@unimore.it), University of Modena and Reggio Emilia (IT)
#' @references Collaud Coen, M., Andrews, E., Bigi, A., Romanens, G., Martucci, G., and Vuilleumier, L.: Effects of the prewhitening method, the time granularity and the time segmentation on the Mann-Kendall trend detection and the associated Sen's slope, Atmos. Meas. Tech., https://doi.org/10.5194/amt-2020-178, 2020.
#' @examples
#'
#' @export


nanprewhite.AR <- function(data, alpha.ak = 95){

    ## first include a modified version of the function
    ## signal::levison to get rid of annoying warning
    durb.levinson <- function (x, p = NULL) 
    {
        fit <- function(acf, p) {
            ref <- numeric(p)
            g <- -acf[2]/acf[1]
            a <- g
            v <- Re((1 - g * Conj(g)) * acf[1])
            ref[1] <- g
            for (t in 2:p) {
                g <- -(acf[t + 1] + a %*% acf[seq(t, 2, by = -1)])/v
                ## original function from the package "signal" modified using as.vector(g) to get rid of
                ## the very annoying warning
                ## "Recycling array of length 1 in array-vector arithmetic"
                a <- c((a + as.vector(g) * Conj(a[seq(t - 1, 1, -1)])), g)
                v <- v * (1 - Re(g * Conj(g)))
                ref[t] <- g
            }
            a <- c(1, a)
            return(list(a = a, v = v, ref = ref))
        }
        if ((!is.null(p) && (p != as.integer(p))) || (p < 2)) 
            stop("p must be integer >= 2.")
        if (is.numeric(x) && is.null(dim(x))) {
            lx <- length(x)
            if (is.null(p) || p >= lx) 
                p <- lx - 1
            r <- fit(x, p)
        }
        else {
            if (is.numeric(x) && !(is.null(dim(x)))) {
                lx <- dim(x)
                if (is.null(p) || p >= lx[1]) 
                    p <- lx[1] - 1
                zr <- apply(x, 2, function(y) fit(y, p))
                zr <- matrix(unlist(zr), nrow = lx[2], byrow = TRUE)
                a <- zr[, 1:(p + 1), drop = FALSE]
                v <- zr[, p + 2]
                ref <- t(zr[, -(1:(p + 2)), drop = FALSE])
                r <- list(a = a, v = v, ref = ref)
            }
            else {
                stop("x must be a numeric vector or matrix.")
            }
        }
        return(r)
    }
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## now the function nanprewhite can start
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    
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
        ## originally it was from signal::levinson
        ## now it is from the function above
        lev  <-  durb.levinson(x1, p)
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
