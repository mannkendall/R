#' Classify data to compute ties
#'
#' compute the number of data considered equivalent and treated as ties
#'
#' @param data this is the vector of data to be analysed
#' @param resolution the measurement resolution, i.e. delta value below which 2 measurements are considered equivalent
#'
#' @return t amount of ties in the data
#'
#' @author Martine Collaud Coen (martine.collaud@meteoswiss.ch), MeteoSwiss (CH) and Alessandro Bigi (abigi@unimore.it), University of Modena and Reggio Emilia (IT)
#' @references Collaud Coen, M., Andrews, E., Bigi, A., Romanens, G., Martucci, G., and Vuilleumier, L.: Effects of the prewhitening method, the time granularity and the time segmentation on the Mann-Kendall trend detection and the associated Sen's slope, Atmos. Meas. Tech., https://doi.org/10.5194/amt-2020-178, 2020.
#' @examples
#'
#' @export

Nb.tie <- function(data, resolution){
    
    if (length(na.omit(data)) > 4){
        m <- min(data, na.rm = T)
        M <- max(data, na.rm = T)
        
        ## determine the domains containing the data

        if ((m < 0 & M < 0) | (m > 0 & M > 0)){
            interval <- abs(M - m)
        } else {
            interval <- M + abs(m)
        }
        
        if (resolution >= interval){
            stop('the given resolution is too large for the considered dataset', call.= TRUE)
        }
        step <- floor(interval / resolution) + 1      
        bins <- seq(from = m, to = m + step * resolution, by = resolution)

        t <- hist( x = data, breaks = bins, plot = FALSE, right = FALSE)$counts
        
    } else {
        t <- NA
    }
    return(t)
}
