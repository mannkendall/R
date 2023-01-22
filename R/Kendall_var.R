#' Estimate of Kendall Variance
#'
#' calculate the variance with ties in the data and ties in time
#'
#' @param data vector of data observations
#' @param t the number of element in each tie
#' @param n number of non missing data for each year, it corresponds to ties in time
#'
#' @return y numeric with the estimate of variance computed as in GAW report No 133 (A. Sirois), p.30 of annex E
#' @author Martine Collaud Coen (martine.collaud@meteoswiss.com), MeteoSwiss (CH) and Alessandro Bigi (abigi@unimore.it), University of Modena and Reggio Emilia (IT)
#' @references Collaud Coen, M., Andrews, E., Bigi, A., Romanens, G., Martucci, G., and Vuilleumier, L.: Effects of the prewhitening method, the time granularity and the time segmentation on the Mann-Kendall trend detection and the associated Sen's slope, Atmos. Meas. Tech., https://doi.org/10.5194/amt-2020-178, 2020.
#'
#' @source  GAW report No 133 (A. Sirois, Environment Canada), eq 4.5 Annex E
#' 
#' @examples
#'
#' @export

Kendall.var <- function(data, t, n){
    
    Lreal <- length(na.omit(data))

    var.t1 <- sum(t * (t - 1) * (2 * t + 5))
    var.t2 <- sum(t * (t - 1) * (t - 2))
    var.t3 <- sum(t * (t - 1))

    var.n1 <- sum(n * (n - 1) * (2 * n + 5))
    var.n2 <- sum(n * (n - 1) * (n - 2))
    var.n3 <- sum(n * (n - 1))

    L <- max(length(data))

    ## message(paste0("L is ", L))
    
    y <- ((( Lreal * (Lreal - 1) * (2 * Lreal + 5)) - var.t1 - var.n1) / 18) +
        (var.t2 * var.n2 / (9 * Lreal * (Lreal - 1) * (Lreal - 2))) +
        (var.t3 * var.n3 / (2 * Lreal * (Lreal - 1)))

    return(y)
}

