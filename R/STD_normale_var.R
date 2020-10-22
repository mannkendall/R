#' Calculate the normal standard variable Z
#' see Gilbert 1987
#' @author Martine Collaud Coen, MeteoSwiss (CH) and alessandro.bigi@unimore.it, University of Modena and Reggio Emilia (IT)
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
