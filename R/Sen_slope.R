#' Sen's slope estimate
#' 
#' It computes the median of the slopes for each interval (xj - xi) / (j - i), j>i (see Gilbert 1987) 
#'     The confidence limits are compute with an interpolation that is important is the number of data is small, such as yearly averages for a 10 year's trend.
#'
#' @param epoch.time: vector of integer times (e.g. days) from the epoch
#' @param data: a vector of observed data
#' @param vari: Kendall variance as calculated by Kendall.var
#' @param alpha.cl: confidence level for the confidence limit (90 or 95%; default: 90%)
#' @return slope: Sen's slope
#' @return LCL: slope lower confidence limit
#' @return UCL: slope upper confidence limit

#' @author Martine Collaud Coen, MeteoSwiss (CH) and alessandro.bigi@unimore.it, University of Modena and Reggio Emilia (IT)
#' @references Collaud Coen, M., Andrews, E., Bigi, A., Romanens, G., Martucci, G., and Vuilleumier, L.: Effects of the prewhitening method, the time granularity and the time segmentation on the Mann-Kendall trend detection and the associated Sen's slope, Atmos. Meas. Tech., https://doi.org/10.5194/amt-2020-178, 2020.
#' @examples
#'
#' @export

Sen.slope <- function(data, epoch.time, vari, alpha.cl = 90){

    data.numeric <- as.numeric(data)
    L <- length(data.numeric)

    dij <- matrix(NA, nrow=L, ncol=L)

    ## compute the slope from all the pairs
    ##(data[i+1:L] - data[i]) / (epoch.time(i+1:L) - epoch.time(i))
    for (i in 1:c(L-1) ){
        dij[ c(i+1):L, i] <- (data.numeric[c(i+1):L] - data.numeric[i]) / (epoch.time[c(i+1):L] - epoch.time[i])
    }

    ## compute the slope from all the pairs
    ## for i=1:L
    ## dij(i+1:L,i) = (data(i+1:L) - data(i)) ./ (time(i+1:L)-time(i));   
    ## end

    ## take the median
    d <- na.omit( as.vector(dij) )
    N <- length(d) 
    slope <- median(d)
    d.ascending <- sort(d)

    ## compute the confidence  limits
    Cconf <- -qnorm( (1-alpha.cl / 100) / 2) * vari ^ 0.5

    M1 <- (N - Cconf) / 2
    M2 <- (N + Cconf) / 2
    if (M1 > 1) {
        ## LCL=d.ascending(round(real(M1)))
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
        ## UCL=d.ascending(round(real(M2)))
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
