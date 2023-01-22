#' Sen's slope estimate
#' 
#' It computes the median of the slopes for each interval (xj - xi) / (j - i), j>i (see Gilbert 1987) 
#'     The confidence limits are compute with an interpolation that is important is the number of data is small, such as yearly averages for a 10 year's trend.
#'
#' @param epoch.time vector of integer times (e.g. days) from the epoch
#' @param data a vector of observed data
#' @param vari Kendall variance as calculated by Kendall.var
#' @param alpha.cl confidence level for the confidence limit (90 or 95; default: 90)
#' @param approx.sol if set to TRUE, Sen's slope is computed by the approximated algorithm by Dillencourt et al. (1992) implemented in the robslopes package. Needed for large datasets.
#' @return slope Sen's slope
#' @return LCL slope lower confidence limit
#' @return UCL slope upper confidence limit
#' 
#' @author Martine Collaud Coen (martine.collaud@meteoswiss.ch), MeteoSwiss (CH) and Alessandro Bigi (abigi@unimore.it), University of Modena and Reggio Emilia (IT)
#' @references Collaud Coen, M., Andrews, E., Bigi, A., Romanens, G., Martucci, G., and Vuilleumier, L.: Effects of the prewhitening method, the time granularity and the time segmentation on the Mann-Kendall trend detection and the associated Sen's slope, Atmos. Meas. Tech., https://doi.org/10.5194/amt-2020-178, 2020.
#' Dillencourt, M. B., Mount, D. M., & Netanyahu, N. S.: A randomized algorithm for slope selection. Inter. J. Comp. Geom. & Applic., https://doi.org/10.1142/S0218195992000020, 1992.
#' @examples
#'
#' @export

sen.slope <- function(data, epoch.time, vari, alpha.cl = 90, approx.sol = FALSE){

    data.numeric <- as.numeric(data)

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    ### Original implementation of the Theil-Sen regression
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

    L <- length(data.numeric)

    ## message("in sen_slope.R: setting dij")
    ## dij <- matrix(NA, nrow=L, ncol=L)

    ## ## compute the slope from all the pairs
    ## message("in sen_slope.R: starting for loop")
    ## for (i in 1:c(L-1) ){
    ##      dij[ c(i+1):L, i] <- (data.numeric[c(i+1):L] - data.numeric[i]) / (epoch.time[c(i+1):L] - epoch.time[i])
    ## }

    ## ## take the median
    ## d <- na.omit( as.vector(dij) )
    ## N <- length(d) 
    ## slope <- median(d)
    ## d.ascending <- sort(d)

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    ### Compute the Theil-Sen regression estimator.
    ### Function by Rand Wilcox, which have been slightly modified
    ### taken through 
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

    temp <- matrix(c(epoch.time, data.numeric), ncol = 2) ## create data matrix
    temp <- temp[-which( is.na(temp[,2]) ), ] ## remove missing values
    x <- temp[, 1]
    y <- temp[, 2]
    
    if(approx.sol == FALSE){

        ord <- order(x)
        xs <- x[ord]
        ys <- y[ord]
        vec1 <- outer(ys, ys, "-")
        vec2 <- outer(xs, xs, "-")
        v1 <- vec1[vec2 > 0]
        v2 <- vec2[vec2 > 0]

        d <- na.omit( as.vector(v1 / v2) )
        N <- length(d) 
        slope <- median(d)
        d.ascending <- sort(d)
        d.ascending <- sort(d)
    }
    
    ## print(paste("MCC: ", N, slope, "RW: ", N1, slope1))

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    ### Compute the Theil-Sen regression estimator.
    ### Use the robslopes package
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    
    if(approx.sol == TRUE) {

        ## message("Computing Sen's slope by numerical approximation")
        
        slope <- robslopes::TheilSen(x, y, alpha = NULL , FALSE)$slope

        ## compute the confidence  limits
        Cconf <- -qnorm( (1 - alpha.cl / 100) / 2) * vari ^ 0.5

        ## number of possible slopes
        N <- length(x) * (length(x) - 1) / 2
        
        M1 <- (N - Cconf) / 2 ## LCL
        M2 <- (N + Cconf) / 2 ## UCL
        
        a1 <- robslopes::TheilSen(x, y, alpha = ceiling(M1)/N, FALSE)$slope
        a2 <- robslopes::TheilSen(x, y, alpha = floor(M1)/N, FALSE)$slope

        LCL <- a2 + (a1 - a2) * (M1 - ifelse(M1 > 0, floor(M1), ceiling(M1)))
        
        a1 <- robslopes::TheilSen(x, y, alpha = ceiling(M2)/N, FALSE)$slope
        a2 <- robslopes::TheilSen(x, y, alpha = floor(M2)/N, FALSE)$slope
 
        UCL <- a2 + (a1 - a2) * (M2 - ifelse(M2 > 0, floor(M2), ceiling(M2)))

        ## message(paste("slope", slope, "M1", M1, "M2", M2, "a11", a11, "a21", a21, "a1", a1, "a2", a2, "LCL", LCL, "UCL", UCL ))
        
        return(list("slope" = slope, "LCL" = LCL, "UCL" = UCL))

    }    

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

    ## compute the confidence  limits
    Cconf <- -qnorm( (1 - alpha.cl / 100) / 2) * vari ^ 0.5

    M1 <- (N - Cconf) / 2
    M2 <- (N + Cconf) / 2

    ## message(paste(M1, N, Cconf, vari))
    
    if (M1 > 1) {
        
        ## interpolation to obtain the best CL
        a1 <- d.ascending[ ceiling(M1) ]
        a2 <- d.ascending[ floor(M1) ]
        LCL <- a2 + (a1 - a2) * (M1 - ifelse(M1 > 0, floor(M1), ceiling(M1)))

    } else if (length(d.ascending == 0)) {
        LCL <- d.ascending[1] ##  the LCL cannot be smaller than the smallest slope
    } else {
        LCL  <- NA
    }
    
    if (M2 > 1){

        ## interpolation to obtain the best CL
        a1 <- d.ascending[ ceiling(M2) ]
        a2 <- d.ascending[ floor(M2) ]
        UCL <- a2 + (a1 - a2) * (M2 - ifelse(M2 > 0, floor(M2), ceiling(M2)))

    } else if (length(d.ascending == 0)) {
        UCL <- d.ascending[ length(d.ascending) ] ##  the UCL cannot be larger than the largest slope
    } else {
        UCL <- NA
    }

    return(list("slope" = slope, "LCL" = LCL, "UCL" = UCL))
}
