#' Applies the 3PW method
#' Estimate the probability of the MK test with the 3PW method. To have the maximal certainty, P is taken as the maximum of P.PW and P.TFPW.Y. Estimate the statistical significance of the MK test as a function of the given confidence level alpha.mk
#'
#' @param P.PW: probability computed from the PW prewhitened dataset
#' @param P.TFPW.Y: probability computed from the TFPW_Y prewhitened dataset
#' @param alpha.mk: confidence level in percentage for the MK test.
#'
#' @return P: probability of the MK test according to the 3PW method
#' @return ss: statistical significance of the MK test (given the alpha.mk)
#' 
#' @author Martine Collaud Coen (martine.collaud@meteoswiss.ch), MeteoSwiss (CH) and Alessandro Bigi (abigi@unimore.it), University of Modena and Reggio Emilia (IT)
#' @references Collaud Coen, M., Andrews, E., Bigi, A., Romanens, G., Martucci, G., and Vuilleumier, L.: Effects of the prewhitening method, the time granularity and the time segmentation on the Mann-Kendall trend detection and the associated Sen's slope, Atmos. Meas. Tech., https://doi.org/10.5194/amt-2020-178, 2020.
#' @examples
#'
#' @export

prob.3PW <- function(P.PW, P.TFPW.Y, alpha.mk) {


    P.alpha <- 1 - alpha.mk / 100
    
    ## compute the probability
    P <- max(P.PW, P.TFPW.Y, na.rm = TRUE)

    ## determine the ss
    if (P.PW <= P.alpha &&  P.TFPW.Y <= P.alpha){
        ss <- alpha.mk
    } else if (P.PW > P.alpha && P.TFPW.Y <= P.alpha) { ## false positive for TFPW.Y @ alpha %
        ss <- -1
    } else if (P.TFPW.Y > P.alpha && P.PW <= P.alpha) { ## false positive for TFPW_Y
        ss <- -2
    } else if (P.TFPW.Y > P.alpha && P.PW > P.alpha) { ## false positive for PW
        ss <-  0
    }

    return(list("P" = P, "ss" = ss))
}
