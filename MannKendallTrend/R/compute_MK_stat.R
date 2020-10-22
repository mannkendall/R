## t.time is array of times expressed as a 6 column data.frame with year, month, day, hour, minute, second
## data is the time series of the variable to be analysed

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
        Prob.MK.n <- read.table('prob_mk_n.csv', sep=",", header = FALSE)
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

    
    out <- Sen.slope( data = data, epoch.time = epoch.time, vari = vari, alpha.cl = alpha.cl)
    result$slope <- out$slope * 365.25
    result$UCL <- out$UCL * 365.25
    result$LCL <- out$LCL * 365.25

    output <- list("result" = result, "S" = S, "vari" = vari, "Z" = Z)
    return(output)
}
