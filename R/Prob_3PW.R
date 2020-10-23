## 1) Estimate the probability of the MK test with the 3PW method.
##    To have the maximal certainty, P is taken as the maximum of P_PW and P_TFPW_Y
## 2) Estimate the statistical significance of the MK test as a function of the
##    given confidence level alpha_MK

## input: 
##       P.PW : probability computed from the PW prewhitened dataset
##       P.TFPW.Y: probability computed from the TFPW_Y prewhitened dataset
##       alpha.mk: level in % for the MK test. 

Prob.3PW <- function(P.PW, P.TFPW.Y, alpha.mk) {


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
