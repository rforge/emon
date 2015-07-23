generate.trend = function(nyears, mu1=0, change, change.type="A", type = c("linear", "incident", "step", "updown"), changeyear, symmetric=F){
  # Generate mean values for different scenarios of change
  # Based on
  # Fryer & Nicholson 1993 IJMS 50:161-168 which had linear change and incident
  # Fryer & Nicholson 1999 IJMS 56:779-790 which had linear, step mid-series, upward then downward
  # although they expressed trends in logs, here not transforming data
  
  # mean values mu_i for years i=1,...,nyears
  # mu1 is value at year 1. Equivalent to fmin in F&N 1999.
  # change is difference between mu1 and largest mu_i, can be negative for a decreasing trend
  # type is type of trend: linear, incident (a single peak), step, updown
  # changeyear defines when incident, step or switch from up to down happens
  
  # for updown: gradient is calculated for upwards change and -1 * same gradient applied for downwards slope
  # special case if symmetric=T and even number of years and changeyear is year T/2
  # then two years in the middle of the series (T/2 and T/2 + 1) both have mu = mu1+k

  type <- match.arg(type) # allows abbreviations for type
   
  if(nyears < 2 | !is.wholenumber(nyears)) stop("Need nyears, the number of years, to be a whole number greater or equal to 2")
  
  if(type!="linear"){
    if(changeyear <1 | changeyear > nyears) stop("Need changeyear to be in range 1 to nyears")
  }

  allowedchange.type = c("M", "A")
  
  if( !(change.type %in% allowedchange.type) )
    stop(paste("change.type must be one of ", allowedchange.type))

  if (change.type=="M") change = mu1*change/100
  
  i = 1:nyears
  
  # using if statements for type, alternatively could use switch but care needed with factors for that
  if(type=="linear"){
    b = change/(nyears-1)
    mu = mu1 + (i-1) * b
  }

  if(type=="incident"){
    mu = rep(mu1, times=nyears)
    mu[changeyear] = mu1 + change
  }
  
  if(type=="step"){
    mu = rep(mu1, times=nyears)
    mu[changeyear:nyears] = mu1 + change
  }
  
  if(type=="updown"){
    if(changeyear == 1 | changeyear == nyears) stop('changeyear cannot be first or last year in series for type="updown". This a linear trend so can use type="linear".')
        
    b = change/(changeyear-1) # slope
    
    if (!(nyears %% 2 == 0 & abs(changeyear - nyears/2)<1e15 & symmetric==T)){    
      mu = mu1 + change - abs(i-changeyear) * b
    } else {
    # special case changeyear in middle of even number of years so trend is symmetrical
    mu = mu1 + (i-1) * b
    mu[(changeyear+1):nyears] = mu[changeyear:1]
    }
  }

  data.frame(i,mu)
}



addnoise = function(meanvalues, reps, distribution="Normal", sd=NA, nbsize=NA,
     randeffect=NA, randeffect.sd=NA){
#***************************************************************************
# Use the generated means as input to a function that adds random variation
# This function is called from within power.trend
#***************************************************************************

nyears = length(meanvalues)
nobs = nyears * reps
meanvalues.rep = rep(meanvalues, rep(reps, nyears))

  
  if(distribution == "Normal")  dat = rnorm(n=nobs, mean = meanvalues.rep, sd = sd)
  if(distribution == "Poisson")  dat = rpois(n=nobs, lambda = meanvalues.rep)  
  if(distribution == "Negbin")  dat = rnbinom(n=nobs, mu = meanvalues.rep, size = nbsize)
  
#  if (randeffect==T) dat = dat + rep(rnorm(n=nyears, mean=0, sd=randeffect.sd), rep(reps, nyears))
  dat
}

power.trend = function(xvalues, reps=1,  meanvalues, distribution="Normal", sd=NA,
  nbsize=NA, method="linear regression", alpha=0.05, nsims=1000, nsims.mk=999,
  randeffect=F, randeffect.sd=NA){
#******************************************************************************
# Power to detect a trend.
#
# Within gam a spline s(xvalues) is used with default settings
# Currently gam only uses Normal errors, could extend this to match distribution in simulation
#*******************************************************************************

  if (missing(meanvalues) || length(meanvalues) == 0L || mode(meanvalues) != "numeric") 
    stop("'meanvalues' must be a non-empty numeric vector")
  if (any(!is.finite(meanvalues))) 
    stop("'meanvalues' contains missing or infinite values")
  
  if(length(distribution)!=1 || mode(distribution) != "character")
    stop('distribution must be a character vector of length one, e.g. "Normal" ')
  
  alloweddist = c("Normal", "Poisson", "Negbin")
  
  if( !(distribution %in%alloweddist) )
    stop(paste("distribution must be one of ", alloweddist))
  
  if(distribution == "Normal" & (missing(sd) || length(sd) == 0L || mode(sd) != "numeric"))
    stop("For distribution Normal, 'sd' must be a non-empty numeric vector")
  
  if(distribution == "Poisson" & any(meanvalues < 0))
    stop("For distribution Poisson 'meanvalues' must be greater or equal to 0")
  
  if(distribution == "Negbin" & (missing(nbsize) || length(nbsize) == 0L || mode(nbsize) != "numeric"))
    stop("For distribution Negbin, 'nbsize' must be a non-empty numeric vector")

  allowedmethod = c("linear regression", "mk", "gam")

  if( !(method %in% allowedmethod) )
    stop(paste("distribution must be one of", allowedmethod))

  if (!is.wholenumber(reps) | reps < 0.9)
      stop("Variable 'reps' must be a positive integer")

  if (mode(randeffect) != "logical")
      stop("The variable 'randeffect' must be logical - either T or F")

if (randeffect & mode(randeffect.sd) != "numeric")
     stop("Random effect sd must be numeric")

if(randeffect & randeffect.sd <=0)
     stop("Random effect sd must be positive")
  
  #if(method == "gam") require("mgcv")

length.xvalues = length(xvalues)
xvalues = rep(xvalues, rep(reps,length.xvalues))

  pvalues = numeric(nsims)

  if(method == "linear regression"){
    for(j in 1:nsims){
    y = addnoise(meanvalues, reps, distribution, sd, nbsize, randeffect, randeffect.sd)
  
    smry = summary(lm(y ~ xvalues))
    pvalues[j] = coef(smry)[2,"Pr(>|t|)"]
    }
  }
  
  if(method == "mk"){
    for(j in 1:nsims){
      y = addnoise(meanvalues, reps, distribution, sd, nbsize, randeffect, randeffect.sd)
      smry = mannkendall(xvalues, y, nsims.mk)
      pvalues[j] = smry$pval
    }
  }

  if(method == "gam"){
    for(j in 1:nsims){
    y = addnoise(meanvalues, reps, distribution, sd, nbsize, randeffect, randeffect.sd)
      
    smry = summary(gam(y~s(xvalues)))
    pvalues[j] = smry$s.table[1,"p-value"]
    }    
  }
  
  # in a test summary gam gave some NaN values
#  hist(pvalues, main="p-values", xlim=c(0,1), breaks=seq(0,1,by=0.05)) 

  # if just report number of NaN values and do calculation without them for the minute
  ispna = is.na(pvalues)
  pvalues.notna = pvalues[!ispna]
  if(sum(ispna)>0)print(paste("For", sum(ispna), "simulations the model summary gave p-value as NA"))
    
  # power
  power = sum(pvalues.notna < alpha) / length(pvalues.notna)
  power
}


mannkendall.stat = function(time, Y) {
#*********************************************************************************
# Calculates Mann-Kendall for a classifying variable (time) and the variable
# of interest (Y). Works even if repeat time points.
# ********************************************************************************
Y.sort = Y[order(time)]
time.sort = time[order(time)]

npoints = length(time)
stat = 0
for (k in 1:(npoints-1)) {
Y.future = Y.sort[time.sort[k] < time.sort]
npos = sum(Y.sort[k] < Y.future)
nneg = sum(Y.sort[k] > Y.future)
stat = stat + npos - nneg
}
stat
}

mannkendall = function(time, Y, nsims.mk=999) {
#****************************************************************************
# Calculates observed mannkendall statistic and also the p-value
# assuming two-sided test (ie alternative hypothesis that trend not equal to
# zero.
# Mann, H.B. (1945), Nonparametric tests against trend, Econometrica, 13,
# 245-259.
# Kendall, M.G. 1975. Rank Correlation Methods, 4th edition, Charles Griffin, 
# London.
#****************************************************************************

# ERROR CHECKS

if(nsims.mk < 10 | !is.wholenumber(nsims.mk)) stop("Number of permutation
   simulations must be a whole number and be at least 10")

if (length(time) != length(Y)) stop("Y and time not of same length")

if (any(!is.finite(time))) 
    stop("time contains missing or infinite values")

if (any(!is.finite(Y))) 
    stop("Y contains missing or infinite values")

pres = !is.na(Y)
Y = Y[pres]; time = time[pres]
pres2 = !is.na(time)
Y = Y[pres2]; time = time[pres2]

observed = mannkendall.stat(time, Y)
npoints = length(time)
theory = rep(0,nsims.mk)

for (j in 1:nsims.mk) {
random = sample(time, npoints, replace=F)
theory[j] = mannkendall.stat(random, Y)
}
bigger = sum(abs(theory) >= abs(observed))
pvalue = (bigger+1)/(nsims.mk+1)

list(pvalue=pvalue, mann=observed)
}



