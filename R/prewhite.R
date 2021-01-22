#' Compute the prewhitened datasets
#'
#' calculate all the necessary prewhitened datasets to asses the statistical significance and to compute the Sen's slope for each of the prewhitening method, including 3PW
#'
#' @param data.ts it is a dataframe with a column named "Time" as POSIXct with timezone UTC and the other columns with observations
#' @param column the number column of the data to be analyzed (must be greater than 2)
#' @param resolution the measurement resolution, i.e. delta value below which 2 measurements are considered equivalent. It is used to compute the number of ties
#' @param alpha.ak statistical significance in percentage for the first lag autocorrelation (default is 95)
#'
#' @return  data.PW a dataframe with 6 columns:
#'      Time (as POSIXct in UTC) = time of the series
#'      PW = PW with the first lag autocorrelation of the data
#'      PW.cor = PW corrected with 1/(1-ak1)
#'      TFPW.WS = PW with the first lag autocorrelation of the data after detrending computed from PW data (see Wang & Swail)
#'      TFPW.Y = method of Yue et al 2002, not 1/1-ak1) correction detrend on original data
#'      VCTFPW = PW with the first lag autocorrelation of the data after detrending + correction of the PW data for the variance (see Wang 2015)
#' 
#' @author Martine Collaud Coen (martine.collaud@meteoswiss.ch), MeteoSwiss (CH) and Alessandro Bigi (abigi@unimore.it), University of Modena and Reggio Emilia (IT)
#' @references Collaud Coen, M., Andrews, E., Bigi, A., Romanens, G., Martucci, G., and Vuilleumier, L.: Effects of the prewhitening method, the time granularity and the time segmentation on the Mann-Kendall trend detection and the associated Sen's slope, Atmos. Meas. Tech., https://doi.org/10.5194/amt-2020-178, 2020.
#' @examples
#'
#' @export

prewhite <- function(data.ts, column, resolution, alpha.ak = 95){

    if (alpha.ak > 100) stop('the confidence limit has to be lower than 100%.');

    ## put the variable to be prewhitened in "data" 
    data <- data.ts[ , column]
    
    ## remove Inf
    data[is.infinite(abs(data))] <- NA

    ## get the vector of days since 1900
    epoch.time <- as.numeric(data.ts$Time) / 3600 / 24

    ## get the dataframe of time split in columns (this is the same as the output of datevec in matlab)
    t.time <- do.call(cbind.data.frame, lapply(c("%Y","%m","%d","%H","%M","%S"), function(x) as.numeric(format(data.ts$Time,x))))
    names(t.time) <- c("year","month","day","hour","minute","seconds")

    ## initialise the loop counter
    nb.loop <- 0

    ## create empty dataframes for postprocessed data
    dataPW <- data.frame("time" = data.ts[,1])
    c <- data.frame(matrix(NA, nrow=1, ncol=0))

    ## computation of the autocorrelation:
    out <- nanprewhite.AR(data = data, alpha.ak = alpha.ak)
    c$PW <- out$ak.lag
    dataARremoved <- as.numeric(out$data.prewhite)
    c$ss <- out$ak.ss
    rm(out)
    
    ##  compute PW-corrected data
    if (length(na.omit(dataARremoved)) > 0 & c$ss == alpha.ak & c$PW >= 0.05) {
        
        dataPW$PW <-  dataARremoved
        dataPW$PW.cor  <- dataARremoved / (1 - c$PW)
                
        ## VCTFPW-corrected data
        ## calculation of the trend slope

        ## compute ties, given the resolution provided
        ties <- Nb.tie(data = dataARremoved / (1 - c$PW), resolution = resolution)
        ## compute the S statistics
        n <- S.test(data = dataARremoved / (1 - c$PW), t.time = t.time)$n
        ## compute the variance for Kendall
        vari <- Kendall.var(data = dataARremoved / (1 - c$PW), t = ties, n = n)

        ## slope of the PW data
        b0.PW <- sen.slope(data = dataARremoved / (1 - c$PW), epoch.time = epoch.time, vari = vari, alpha.cl = 90)$slope ##
        ## compute ties and the slope for the original data
        t  <- Nb.tie(data = data, resolution = resolution)
        n <- S.test(data = data, t.time = t.time)$n
        vari <- Kendall.var(data = data, t = ties, n = n)
        b0.or <- sen.slope(data = data, epoch.time = epoch.time, vari = vari)$slope

        ## remove the trend
        dataDetrend.PW <- data - b0.PW * c(epoch.time - epoch.time[1])
        dataDetrend.or <- data - b0.or * c(epoch.time - epoch.time[1])

        ## compute the autocorrelation of the de-trended time series:       
        ## [c.VCTFPW, dataARremovedor, c.ssVC] = nanprewhite_ARok(dataDetrendor,'alpha_ak',alpha_ak);
        out <-  nanprewhite.AR(data = dataDetrend.or, alpha.ak = alpha.ak)
        c$VCTFPW <- out$ak.lag
        dataARremoved.or <- as.numeric(out$data.prewhite)
        c$ssVC <- out$ak.ss
        rm(out)

        ## [akPW, dataARremovedPW, ssPW]= nanprewhite.AR(dataDetrend.PW,'alpha_ak',alpha_ak);
        out <- nanprewhite.AR(data = dataDetrend.PW, alpha.ak = alpha.ak)
        ak.PW <- out$ak.lag
        dataARremoved.PW <- as.numeric(out$data.prewhite)
        ss.PW <- out$ak.ss
        rm(out)
        
        ## computation of TFPW correction Yue et al., 2002
        ## blended data
        if (sum( na.omit(dataARremoved.or) > 0)) {
            dataPW$TFPW.Y <-  dataARremoved.or + b0.or * (epoch.time - epoch.time[1])
        } else {
            dataPW$TFPW.Y <- data
        }
        
        ## computation of TFPW correction Wang and Swail
    
        if (abs(ak.PW) >= 0.05 & ss.PW == 95){

            ## change so that while can be used with the same variable:
            ## ak is new c and c1 is last c
            ## ak = c.VCTFPW
            
            c1 <- c$PW
        
            dataARremoved.PW  <-  c(data[1], (data[-1] - ak.PW * data[-length(data)]) / (1 - ak.PW))
        
            t <- Nb.tie(data = dataARremoved.PW, resolution = resolution)

            n <- S.test(data = dataARremoved.PW, t.time = t.time)$n
            
            vari <- Kendall.var(dataARremoved.PW, t = ties, n = n)
            
            b1.PW <- sen.slope(data = dataARremoved.PW, epoch.time = epoch.time, vari = vari)$slope
            
            ## remove the trend

            while  (abs(ak.PW - c1) > 0.0001 & abs(b1.PW - b0.PW) > 0.0001) {
                if (ak.PW >= 0.05 & ss.PW == 95) {
                    nb.loop <- nb.loop + 1
                    dataDetrend.PW <- data - b1.PW * (epoch.time - epoch.time[1])

                    c1 <- ak.PW
                    b0.PW <- b1.PW

                    out <-  nanprewhite.AR(data = dataDetrend.PW, alpha.ak = alpha.ak)
                    ak.PW <- out$ak.lag
                    dataARremoved.2PW <- as.numeric(out$data.prewhite)
                    ss.PW <- out$ak.ss
                    rm(out)
                    
                    if (ak.PW > 0 & ss.PW == 95) {

                        dataARremoved.2PW <- c( data[1], c(data[-1] - ak.PW * data[-length(data)]) / (1 - ak.PW) )

                        t <- Nb.tie(data = dataARremoved.2PW, resolution = resolution)
                        n <- S.test(data = dataARremoved.2PW, t.time = t.time)$n
                        vari <- Kendall.var(data = dataARremoved.2PW, t = t, n = n)
                        b1.PW <- sen.slope(data = dataARremoved.2PW, epoch.time = epoch.time, vari = vari)$slope
                        dataARremoved.PW <- dataARremoved.2PW
                        if (nb.loop > 10) break
                        
                    } else {
                        break
                    }
                }
            }
        } else {
            b1.PW <- b0.PW
            ak.PW <- c$PW
        }

        ## blended data
        if (sum ( na.omit(dataARremoved.PW) ) > 0) {
            dataPW$TFPW.WS <- dataARremoved.PW
            c$TFPW.WS <- ak.PW
        } else {
            dataPW$TFPW.WS <- data
        }
    
        ## correction VCTFPW
        ## correction of the variance
        var.data <- var(data, na.rm = TRUE)
        var.data.TFPW <- var(dataARremoved.or, na.rm = TRUE)
        dataARremoved.var <- dataARremoved.or * (var.data / var.data.TFPW)

        ## modify slope estimator (correction of the slope for the autocorrelation)
        if (c$VCTFPW >= 0){
            bVC <- b0.or / sqrt( (1 + c$VCTFPW) / (1 - c$VCTFPW) )
        } else if(c$VCTFPW < 0){
            bVC <- b0.or
        }
        ## add the trend again
        dataPW$VCTFPW <-  dataARremoved.var + bVC * (epoch.time - epoch.time[1])
        
    } else {
        ## no s.s. autocorrelation
        dataPW$PW <-  data
        dataPW$PW.cor <-  data
        dataPW$TFPW.Y <-  data
        dataPW$TFPW.WS <-  data
        dataPW$VCTFPW <-  data
    }
    
    return(dataPW)
}
