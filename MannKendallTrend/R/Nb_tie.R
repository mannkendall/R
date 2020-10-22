##  compute the number of data considered as equivalent and treated as ties
##  input: "data" is the vector of data to be analysed
##         "resolution" is the measurement resolution, i.e. the interval between which 2 measures are considered as equals
##         "data_tot" is necessary so that all columns have the same lenght. It is used to compute ties.
##
##  output: "vector" with the number of data per tie 

## tested Oct 12th ok

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

        ## Then compute the number of elements in each bin.
        ## return np.histogram(data, bins=bins)[0]

        ## python-based
        ## nbins  <- as.integer( floor((M - m) / resolution) + 1)
        ## nbins = int((np.nanmax(data)-np.nanmin(data))//resolution + 1)
        ## bins2 <- seq(from = m, to = m + step * resolution, length.out = nbins + 1)
        ## bins = np.linspace(np.nanmin(data), np.nanmin(data) + nbins * resolution, num=nbins + 1)

        ## matlab-based
        ## this is called step, but actually it's the number of bins
        step <- floor(interval / resolution) + 1      
        bins <- seq(from = m, to = m + step * resolution, by = resolution)

        t <- hist( x = data, breaks = bins, plot = FALSE, right = FALSE)$counts
        ## t <- hist( x = data, breaks = bins, plot = FALSE, right = FALSE)
        ## t2 <- hist( x = data, breaks = bins2, plot = FALSE, right = FALSE)$counts
        
    } else {
        t <- NA
    }
    return(t)
}
