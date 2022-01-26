#' compute Mann-Kendall statistics
#'
#' This function computes Mann Kendall statistics for a given time series
#'
#' @param t.time is array of times expressed as a 6 column data.frame with year, month, day, hour, minute, second
#' @param data is the time series of the variable to be analysed
#' @param resolution the measurement resolution, i.e. delta value below which 2 measurements are considered equivalent. It is used to compute the number of ties
#' @param alpha.mk confidence limit for Mk test in percentage. Default value is 95
#' @param alpha.cl confidence limit for the confidence limits of the Sen's slope in percentage. Default value is set to 90
#'
#' @return a list named `output` containing a data.frame (`result`) and 3 numerics (`S`, `vari`, `Z`). `result` is a dataframe with 3 columns: `slope`, the estimate of the Sen's slope (as \% year^-1), `UCL` and `LCL` respectively the upper and the lower confidence limits of the slope (as \% year^-1). `S` is the value of the S statiscs, `vari` is the Kendall variance, `Z` is the resulting standard normal variable
#' 
#' @author Martine Collaud Coen (martine.collaud@meteoswiss.ch), MeteoSwiss (CH) and Alessandro Bigi (abigi@unimore.it), University of Modena and Reggio Emilia (IT)
#' @references Collaud Coen, M., Andrews, E., Bigi, A., Romanens, G., Martucci, G., and Vuilleumier, L.: Effects of the prewhitening method, the time granularity and the time segmentation on the Mann-Kendall trend detection and the associated Sen's slope, Atmos. Meas. Tech., https://doi.org/10.5194/amt-2020-178, 2020.
#' @examples
#'
#' @export
#' @importFrom magrittr %>%

compute.MK.stat <- function(data, t.time, resolution, alpha.mk = 95, alpha.cl = 90) {

    t <- Nb.tie(data = data, resolution = resolution)
    out <- S.test(data = data, t.time = t.time)
    S <- out$S
    n <- out$n
    result <- list()

    vari <- Kendall.var(data = data, t = t, n = n)
    Z <- STD.normale.var(data = S, var.data = vari)
    if ( sum(data, na.rm = TRUE) > 10 ) {
        result$P <- 2 * (1 - pnorm( abs(Z), 0, 1))
    } else {
        ## Prob.MK.n <- read.table('prob_mk_n.csv', sep=",", header = FALSE)
        result$P <- Prob.MK.n[ abs(S) + 1, sum(data, na.rm = TRUE)]
    }
    
    ## determine the ss
    if (result$P <= c(1 - alpha.mk/100)){
        result$ss= alpha.mk
    } else {
        result$ss <- 0
    }

    apply(t.time, 1, function(x) paste0(x[1:6], collapse="-")) %>%
        unlist() %>%
        strptime(., format= c("%Y-%m-%d-%H-%M-%S"), tz = "UTC") %>%
        as.POSIXct(., tz="UTC") %>%
        as.numeric(.) / 3600 /24 -> epoch.time

    
    out <- sen.slope( data = data, epoch.time = epoch.time, vari = vari, alpha.cl = alpha.cl)
    result$slope <- out$slope * 365.25
    result$UCL <- out$UCL * 365.25
    result$LCL <- out$LCL * 365.25

    output <- list("result" = result, "S" = S, "vari" = vari, "Z" = Z)
    return(output)
}
