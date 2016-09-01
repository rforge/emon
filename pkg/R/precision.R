
# function to check for whole numbers, taken from is.integer help file
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

n.min = function(N, sigma, d, prop) {
  #*****************************************************************************
  # Used with optimize() because the t-value is also a function of N
  #*****************************************************************************
  f = abs(N - 4*(qt(prop,(N-1)) * sigma / d)^2)
  f
}


precision = function(d, n, pars, method="sample size", alpha=0.05,
                     minint=1, maxint=500)  {
  #***********************************************************************
  # Calculates sample size needed for precision d (width of confidence
  # interval for the mean) of a single sample.
  #***********************************************************************
  #
  #### Parameter testing ####
  
  allowedmethods = c("sample size", "width")
  if( !(method %in% allowedmethods) )
    stop("Tests must be of type character and be one of sample size
         or width")
  
  if (alpha<=0 | alpha>=1)
    stop("alpha must be a positive real number less than 1")
  
  if (pars<=0)
    stop("pars must be a positive real number")
  
  if (method == "sample size") {
    length.d = length(d)
    for (k in 1:length.d) {
      if (d[k]<0) stop("d must be a positive real number")
    }}
  
  if (method == "width") {
    length.n = length(n)
    
    for (k in 1:length.n) {
      if ( !is.wholenumber(n[k]) | n[k]<1.5 ) 
        stop("n must be an integer of at least 2")
    }}
  
  #### End of parameter testing #### 
  
  if (method == "sample size") {
    n = vector("numeric", length.d)
    
    sigma1 = pars[1]
    prop = 1 - alpha/2
    
    for (k in 1:length.d) {
      ssize = optimize(n.min, interval=c(minint, maxint), sigma=sigma1, d=d[k],
                       prop=prop, maximum=F)$min
      n[k] = ceiling(ssize)
      if (n[k]>=maxint) warning("maxint reached", call. = FALSE)
      if (n[k]<=(minint+1)) warning("minint reached", call. = FALSE)
    }
  }
  
  if (method == "width") {
    width = vector("numeric", length.n)
    sigma = pars[1]
    for (k in 1:length.n) {
      prop = 1 - (alpha/2)
      width[k] = 2 * sigma * qt(prop, n[k]-1) / sqrt(n[k])
    }
    d = width
  }
  list(n = n, d = d)
}

