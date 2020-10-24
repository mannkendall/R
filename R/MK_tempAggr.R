#' MK test and Sen slope calculator
#' 
#' MK test and the Sen slope on the given time granularity. The function implements three prewhitening (PW) methods: PW (Yue et al., 2002) and TFPW.Y (trend free PW, Wang and Swail, 2001) to compute the statistical significance and VCTFPW (Wang, W., Chen, Y., Becker, S., & Liu, B. (2015). Variance Correction Prewhitening Method for Trend Detection in Autocorrelated Data. J. Hydrol. Eng., 20(12), 4015033-1-04015033�10. https://doi.org/10.1061/(ASCE)HE.1943-5584.0001234.) to compute the Sen's slope. Only the statistically significant (ss) autocorrelation are taken into account for the prewhitening. The ss of the trends is taken at 95% confidence limit. The upper and lower confidence limits are given by the 90% of the all intervals differences distribution. The significance level is given by the MK test and has therefore no direct relation to the confidences limits. If seasonal Mann-Kendall is applied, the yearly trend is assigned only if the results of the seasonal test are homogeneous. 
#' 
#' @references: WMO-GAW publication N. 133, annex E, p. 26, MULTMK/PARTMK by C. Libiseller and the book of Gilbert 1998
#' @references: Collaud Coen, M., Andrews, E., Bigi, A., Romanens, G., Martucci, G., and Vuilleumier, L.: Effects of the prewhitening method, the time granularity and the time segmentation on the Mann-Kendall trend detection and the associated Sen's slope, Atmos. Meas. Tech., https://doi.org/10.5194/amt-2020-178, 2020.
#' 
#' @param data: a data.frame with the first column being a POSIXct object with tz = "UTC" and the second column the variable to be analysed
#' @param PW.method: the PW method used: e.g. one among: 3PW, PW, TFPW.Y, TFPW.WS and VCTFPW. The default is 3PW.
#' @param resolution: It is taken into account to determine the number of ties; a good guess it is the resolution of the instrument or a little bit higher. This parameters can be determinant for the results but not very sensitive
#' @param alpha.mk: confidence limit for Mk test in percentage. Default value is 95
#' @param alpha.cl: confidence limit for the confidence limits of the Sen's slope in percentage. Default value is 90 
#' @param alpha.xhomo: confidence limit for the homogeneity between seasons in percentage. Default value is 90 
#' @param alpha.ak: confidence limit for the first lag autocorrelation in percentage. Default value is 95
#' @param seasonal: set to TRUE if the analysis needs to be performed over user-defined seasons (default is FALSE)
#' @param seasons: a vector of the same lenght of the number of records in data.tempAgg used to split the data into seasons. It is used only if seasonal = TRUE.  
#' @return P: probability for the statistical significance. If 3PW is applied, P = max(P.PW, P.TFPW.Y)
#' @return ss: statistical significance: alpha % if the test is ss at the alpha confidence level. Default = 95%. 0 if the test is not ss at the alpha confidence level; -1 if the test is a TFPW.Y false positive at alpha confidence level; -2 if the test is a PW false positive at alpha confidence level
#' @return slope: Sen's slope in units/y
#' @return UCL: upper confidence level in units/y
#' @return LCL: lower confidence level in units/
#' 
#' @author Martine Collaud Coen, MeteoSwiss (CH) and alessandro.bigi@unimore.it, University of Modena and Reggio Emilia (IT)
#' @references Collaud Coen, M., Andrews, E., Bigi, A., Romanens, G., Martucci, G., and Vuilleumier, L.: Effects of the prewhitening method, the time granularity and the time segmentation on the Mann-Kendall trend detection and the associated Sen's slope, Atmos. Meas. Tech., https://doi.org/10.5194/amt-2020-178, 2020.
#' @examples
#'
#' @export

MK.tempAggr <- function(data, PW.method = "3PW", resolution, alpha.mk = 95, alpha.cl = 90, alpha.xhomo = 90, alpha.ak = 95, seasonal = FALSE, seasons = NULL){

    data <- data[order(data$Time), ]
   
    ## check the arguments
    if( sum(PW.method %in% c("3PW", "PW", "TFPW_Y", "TFPW_WS", "VCTFPW")) == 0) stop("Invalid argument for the PW method")

    if (seasonal == FALSE){

        ## prepare result data
        result <- data.frame(matrix(NA, nrow = 1, ncol = 5))
        names(result)  <- c("slope","UCL","LCL","P","ss")

        ## dataframe with 6 columns for year, month, day, hour, minute, second 
        t.time <- lapply(c("%Y","%m","%d","%H","%M", "%S"), function(x) as.numeric(format(data[, 1], x))) %>%
            as.data.frame(.)
        
        ## compute all the prewhitened time series
        ## input of prewhite is data.frame with a first posix column
        ## dataPW is a dataframe containing the 3 prewhitened series
        ## data.ak_y is the first lag autocorrelation coefficient for the complete time series

        message("Prewithening the timeseries...")

        dataPW <- prewhite(data.ts = data, column = 2, resolution = resolution, alpha.ak = alpha.ak)

        pw.met.col <- grep(PW.method, names(dataPW), fixed = TRUE)

        if (is.na(pw.met.col) == FALSE && PW.method != "3PW") {
            d.comp <- data[ , pw.met.col]
            
            ## compute the  Mann-Kendall test without temporal aggregation
            result <- compute.MK.stat(data = d.comp, t.time = t.time, resolution = resolution, alpha.mk = alpha.mk, alpha.cl = alpha.cl)$result
        }

        if (PW.method == '3PW')  {

            result.PW <- compute.MK.stat(data = dataPW$PW, t.time = t.time, resolution = resolution, alpha.mk = alpha.mk, alpha.cl = alpha.cl)$result

            result.TFPW.Y <- compute.MK.stat(data = dataPW$TFPW.Y, t.time = t.time, resolution = resolution, alpha.mk = alpha.mk, alpha.cl = alpha.cl)$result

            result.TFPW.WS <- compute.MK.stat(data = dataPW$TFPW.WS, t.time = t.time, resolution = resolution, alpha.mk = alpha.mk, alpha.cl = alpha.cl)$result

            result.VCTFPW <- compute.MK.stat(data = dataPW$VCTFPW, t.time = t.time, resolution = resolution, alpha.mk = alpha.mk, alpha.cl = alpha.cl)$result

            ## determine the P and ss
            out <- prob.3PW(P.PW = result.PW$P, P.TFPW.Y = result.TFPW.Y$P, alpha.mk = alpha.mk)
            
            result$P <- out$P
            result$ss <- out$ss
            result$slope <- result.VCTFPW$slope
            result$UCL <- result.VCTFPW$UCL
            result$LCL <- result.VCTFPW$LCL
            
        }
    }

    if (seasonal == TRUE){
 
        message("Prewithening the full timeseries...")
        
        ## First prewhite the full timeseries
        dataPW <- prewhite(data.ts = data, column = 2, resolution = resolution, alpha.ak = alpha.ak)

        ## split the time series across user-defined seasons
        data.season <- split(x = dataPW, f = seasons)

        ## prepare result dataframe
        result <-  data.frame(matrix(NA, nrow = length(data.season) + 1, ncol = 5))
        colnames(result)  <- c("slope","UCL","LCL","P","ss")

        ## results needed for non-3PW method
        Z <- rep(NA, length(data.season))
        Stot <- NA
        vari.tot <- NA

        ## results needed for 3PW method
        S.PW <- NA
        vari.PW <- NA
        S.TFPW.Y <- NA
        vari.TFPW.Y <- NA
        Z.VCTFPW <- rep(NA, length(data.season))
        
        ## loop through the elements of data.season (i.e. each season)
        for (m in 1:length(data.season)) {

            message("Analysing season ", m, "...")

            data.mois <- data.season[[m]]
            
            ## get the matrix of times as year, month, day, hour, minute
            t.time.mois <- lapply(c("%Y","%m","%d","%H","%M", "%S"), function(x) as.numeric(format(data.mois$time, x))) %>%
                as.data.frame(.)

            if ( sum( na.omit(data.mois$PW) ) > 1) { ## un peu trop peu ??

                ## get the correct column in data.mois depending on PW.method
                pw.met.col <- grep(PW.method, names(data.mois), fixed = TRUE)

                if ( is.na(pw.met.col) == FALSE && PW.method != "3PW") {

                    ## get the seasonal data
                    d.comp <- data.mois[ , pw.met.col]
                   
                    ## compute the  Mann-Kendall test using the PW.method
                    out <- compute.MK.stat(data = d.comp, t.time = t.time.mois, resolution = resolution, alpha.mk = alpha.mk, alpha.cl = alpha.cl)
                    
                    result[m, c("slope","UCL","LCL")] <- out$result[c("slope","UCL","LCL")]

                    ## update Z, total S statistics, total variance
                    Z[m] <- out$Z
                    Stot <- sum(Stot, out$S, na.omit = TRUE)
                    vari.tot <- sum(vari.tot, out$vari, na.omit = TRUE)
                }

                if (PW.method == "3PW") {

                    out.PW <- compute.MK.stat(data = data.mois$PW, t.time = t.time.mois, resolution = resolution, alpha.mk = alpha.mk, alpha.cl = alpha.cl)

                    out.TFPW.Y <- compute.MK.stat(data = data.mois$TFPW.Y, t.time = t.time.mois, resolution = resolution, alpha.mk = alpha.mk, alpha.cl = alpha.cl)
                    
                    out.TFPW.WS <- compute.MK.stat(data = data.mois$TFPW.WS, t.time = t.time.mois, resolution = resolution, alpha.mk = alpha.mk, alpha.cl = alpha.cl)
                    
                    out.VCTFPW <- compute.MK.stat(data = data.mois$VCTFPW, t.time = t.time.mois, resolution = resolution, alpha.mk = alpha.mk, alpha.cl = alpha.cl)

                    out <- prob.3PW("P.PW" = out.PW$result$P, "P.TFPW.Y" = out.TFPW.Y$result$P, "alpha.mk" = alpha.mk)
                    
                    ## update result table for the m-th season
                    result[m, c("P","ss")] <- out[c("P","ss")]
                    result[m, c("slope", "UCL", "LCL")]  <- out.VCTFPW$result[c("slope", "UCL", "LCL")]

                    S.PW <- sum(S.PW, out.PW$S, na.rm = TRUE)
                    vari.PW <- sum(vari.PW, out.PW$vari, na.rm = TRUE)
                    
                    S.TFPW.Y <- sum(S.TFPW.Y, out.TFPW.Y$S, na.rm = TRUE)
                    vari.TFPW.Y <- sum(vari.TFPW.Y, out.TFPW.Y$vari, na.rm = TRUE)
                    Z.VCTFPW[m] <- out.VCTFPW$Z

                }
                
            } else {
                ## there is no sufficient data for this season
                result[m, c("P", "ss", "slope", "UCL", "LCL")] <- NA
            }
        }

        
        ## compute for the whole period, that is the whole year
        message("Testing for seasonal homogeneity...")

        if (PW.method %in% c('PW','TFPW_Y','TFPW_WS','VCTFPW')) {

            Ztot <- STD.normale.var( Stot, vari.tot )

            if (sum( dataPW$PW, na.rm = TRUE) > 10) {
                result[m+1, "P"] <- 2 * (1 - pnorm( abs(Ztot), 0, 1))
            } else {
                ## Prob.MK.n <- read.table('prob_mk_n.csv', sep=",", header = FALSE)
                result[m + 1, "P"] <- Prob.MK.n[ abs(Stot) + 1, sum( dataPW$PW, na.rm = TRUE)]
            }

            if (result[m + 1, "P"] <= 1 - alpha.mk / 100) {
                result[m + 1, "ss"] <-  alpha.mk
            } else {
                result[m + 1, "ss"] <- 0
            }

            ## compute Xi-square to test the homogeneity between months.
            ## TO DO: change Xhomo as a function of alpha_Xhomo and the number of
            ## the seasons
            Xhomo <- sum ( na.omit(Z)^2 ) - 12 * (mean(Z, na.rm = TRUE))^2
        }

        if (PW.method == "3PW") {

            ## compute the statistical significance for PW
            Ztot.PW <- STD.normale.var( S.PW, vari.PW )

            if (sum( dataPW$PW, na.rm = TRUE) > 10) {
                Ptot.PW <- 2 * (1 - pnorm( abs(Ztot.PW), 0, 1))
            } else {
                ## Prob.MK.n <- read.table('prob_mk_n.csv', sep=",", header = FALSE)
                Ptot.PW <- Prob.MK.n[ abs(Stot.PW) + 1, sum( dataPW$PW, na.rm = TRUE) ]
            }

            ## compute the statistical significance for TFPW_Y
            Ztot.TFPW.Y <- STD.normale.var( S.TFPW.Y, vari.TFPW.Y)
            if ( sum(dataPW$TFPW.Y, na.rm = TRUE) > 10){
                Ptot.TFPW.Y <- 2 * (1 - pnorm( abs(Ztot.TFPW.Y), 0, 1))
            } else {
                ## Prob.MK.n <- read.table('prob_mk_n.csv', sep=",", header = FALSE)
                Ptot.TFPW.Y <- Prob.MK.n[ abs(S.TFPW.Y) + 1, sum( dataPW$TFPW.Y, na.rm = TRUE)] 
            }
            
            ## determine the ss
            out <- prob.3PW(P.PW = Ptot.PW, P.TFPW.Y = Ptot.TFPW.Y, alpha.mk = alpha.mk)
            result[m+1, c("P", "ss")] <- c(out$P, out$ss)
            
            ## compute xi-square to test the homogeneity between months. Since the slope
            ## is computeated from VCTFPW, the homogeneity is also computeated from VCTFPW
            Xhomo <- sum(Z.VCTFPW^2, na.rm = TRUE) - 12 * (mean(Z.VCTFPW, na.rm = TRUE))^2
        }
        
        ## write the yearly slope and CL
        ## Xhomo has a chi-squared distributions with n-1 and 1 degree of
        ## freedom. Seasonal trends are homogeneous is Xhomo is smaller
        ## than the threshold defined by the degree of freedom and the
        ## confidence level alpha_Xhomo.
        ## change condition: yearly slope not given if the seasons are not
        ## homogeneous
        
        if (Xhomo <= qchisq("p" = 1 - alpha.xhomo / 100, "df" = length(data.season) - 1)) {
            result$slope[m + 1] <- median( result$slope[1:m], na.rm = TRUE )
            result$UCL[m + 1] <- median( result$UCL[1:m], na.rm = TRUE )
            result$LCL[m + 1] <- median( result$LCL[1:m], na.rm = TRUE )
        } else {
            warning('the trends for the temporal aggregation are not homogeneous')
            result$slope[m + 1] <- NA
            result$UCL[m + 1] <- NA
            result$LCL[m + 1] <- NA
        }
    }
    return(result)
}
