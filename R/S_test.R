#' Compute the S statistic for the MK test
#'
#' 
#' @param data it is the vector of observations
#' @param t.time it is a matrix with the year, month, date vector as individual variables in each column. The number of rows is the number of obsverations.
#'
#' @return a list containing S and n. S is the double sum on the sign of the difference between data pairs (Si) and n is the number of valid data in each year of the time series
#' 
#' @author Martine Collaud Coen, MeteoSwiss (CH) and alessandro.bigi@unimore.it, University of Modena and Reggio Emilia (IT)
#' @references Collaud Coen, M., Andrews, E., Bigi, A., Romanens, G., Martucci, G., and Vuilleumier, L.: Effects of the prewhitening method, the time granularity and the time segmentation on the Mann-Kendall trend detection and the associated Sen's slope, Atmos. Meas. Tech., https://doi.org/10.5194/amt-2020-178, 2020.
#' @examples
#'
#' @export

S.test <- function(data, t.time){
    
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
        n[i - an.debut + 1] <- length( na.omit(data[t.time[,1] == i]) )

        if (i < an.end){

            for (j in c(i + 1) : an.end) {

                nj <- length(data[t.time[,1]==j])
                sj <- rep(NA, nj)
            
                dataj <- data[t.time[,1]==j]

                for (k in 1:nj){
                    sj[k] <- sum( sign( dataj[k] - data[t.time[,1] == i]), na.rm=TRUE)
                }
                
                Sij[j - an.debut + 1, i - an.debut + 1] <- sum(sj, na.rm=TRUE)
            }
        }
    }
    
    S <- sum(Sij, na.rm = TRUE)   
    return(list("n"=n, "S"=S))
}
