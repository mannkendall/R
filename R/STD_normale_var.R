#' Compute the normal standard variable Z
#'
#' Compute the standard normal variable Z
#'
#' @param data variable to be analysed
#' @param var.data variance of the data
#' @return Z the normal standard variable
#' 
#' @references Gilbert, R. O.: Statistical Methods for Environmental Pollution Monitoring, Van Nostrand Reinhold Company, New York, 1987
#' @author Martine Collaud Coen (martine.collaud@meteoswiss.ch), MeteoSwiss (CH) and Alessandro Bigi (abigi@unimore.it), University of Modena and Reggio Emilia (IT)
#' 
#' @examples
#'
#' @export

STD.normale.var <- function(data, var.data){
    if (data == 0) Z  <-  0
    
    if (data > 0){
        Z <- (data - 1) / ((var.data) ^ 0.5)
    } else {
        Z <- (data + 1) / ((var.data) ^ 0.5)
    }
    
return(Z)
}
