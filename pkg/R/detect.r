
# function to check for whole numbers, taken from is.integer help file
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


fS.detect = function(X, pdet) {
#*****************************************************************************
# Used to find sample size for Square Lattice design, with optimize()
#*****************************************************************************
f =  abs(pdet - (X^2 * (pi - 4*acos(1/(2*X))) + sqrt(4*X^2 - 1)))
f }


fT.detect = function(X, pdet) {
#*****************************************************************************
# Used to find sample size for Triangular Lattice design, with optimize()
#*****************************************************************************
bound1 = 2^(-0.5) * 3^(-0.25)
Z = acos(bound1 / X)
f = abs(pdet - (pi*X^2 - 6*Z*X^2 + sqrt(6*sqrt(3)*X^2 - 3)))
f }


detect = function(method, statistic, area=NA, radius=NA, pdetect=NA, ssize=NA) {
#*****************************************************************************
# Updated by JB 16/9/14
# Updayed by JB 23/2/15 - RANDOM DESIGN FOR FIXED N
# Calculates various design statistics for Random, Square and
# Triangular spatial designs
#
# Method    = 'R' (random)
#           = 'S' (square)
#           = 'T' (triangular)
# Statistic = 'P' (probability detection)
#           = 'N' (sample size)
#           = 'R' (patch radius)
#
# Output    $prob  = probability patch detected
#           $ssize = sample size
#           $rad   = patch radius
#           $sep   = separation distance (for square and triangular designs)
#****************************************************************************
#***********************************************
# Error Checks on parameters
#***********************************************
allowedmethods = c("R", "S", "T")
if( !(method %in% allowedmethods) )
#     stop(paste("Sampling pattern must be one of", allowedmethods))
      stop("Sampling pattern must be of type character and be one of R, S or T")

allowedstatistics = c("P", "N", "R")
if( !(statistic %in% allowedstatistics) )
      stop("Statistic must be of type character and one of P, N or R")

#************************************************************************
# Check that sufficient parameters have been input
#************************************************************************
if (statistic=="P") {
if (is.na(area)==T | is.na(radius)==T | is.na(ssize)==T) {
stop("area, radius and ssize arguments needed")
}}

if (statistic=="N") {
if (is.na(area)==T | is.na(radius)==T | is.na(pdetect)==T) {
stop("area, radius and pdetect arguments needed")
}}

if (statistic=="R") {
if (is.na(area)==T | is.na(ssize)==T | is.na(pdetect)==T) {
stop("area, ssize and pdetect arguments needed")
}}

#************************************************************************

if (statistic != "N") {
if ( !is.wholenumber(ssize) | ssize<0.5 )
      stop("Sample size must be a positive integer")
}

if (area <=0) stop("area must be positive")

if (statistic != "R") { 
if (radius <=0) stop("radius must be positive")
}

if (statistic != "P" & method == "R") {
if (pdetect <=0 | pdetect >=1)
       stop("pdetect must be greater than 0 and less than 1")
}

if (statistic != "P" & (method == "S" | method == "T")) {
if (pdetect <=0 | pdetect >1)
       stop("pdetect must be greater than 0 and not greater than 1")
}

#****************************

if (statistic == 'N') {
patcharea = pi * radius^2
if ( patcharea >= area )
      stop("Patch area must be smaller than sampling area")
}

#**************************
# Random design
#**************************
if (method == 'R') {
if (statistic == 'P') {
density = area / ssize
d = sqrt(density)
R = radius / d
pdetect = 1 - (1-(pi*radius*radius / area))^ssize
}
if (statistic == 'N') {
ssize = log(1-pdetect) / log(1-(pi*radius*radius / area))
}
if (statistic == 'R') {
radius = sqrt(area*(1 - (1-pdetect)^(1/ssize))/pi)
}
distance = NA
}
#**************************
# Square lattice design
#**************************
if (method == 'S') {
if (statistic == 'P') {
density = area / ssize
d = sqrt(density)
R = radius / d
if (R <= 0.5) pdetect = pi * R^2
if ((R > 0.5) & (R <= sqrt(0.5))) pdetect = R^2 * (pi - 4*acos(1/(2*R))) + sqrt(4*R^2 - 1)
if (R > sqrt(0.5)) pdetect = 1
}
if (statistic == 'N') {
pdetect.1 = pi * 0.5^2
if (pdetect <= pdetect.1) R = sqrt(pdetect / pi)
if ((pdetect > pdetect.1) & (pdetect <= 1)) R = optimize(fS.detect, interval=c(0.5, sqrt(0.5)), pdet = pdetect, maximum=F)$min
if (pdetect == 1) R = sqrt(0.5)
d = radius / R
ssize = area / d^2
}
if (statistic == 'R') {
pdetect.1 = pi * 0.5^2
if (pdetect <= pdetect.1) R = sqrt(pdetect / pi)
if ((pdetect > pdetect.1) & (pdetect <= 1)) R = optimize(fS.detect, interval=c(0.5, sqrt(0.5)), pdet = pdetect, maximum=F)$min
if (pdetect == 1) R = sqrt(0.5)
d = sqrt(area / ssize)
radius = R * d
}
distance = d
}
#**************************
# Triangular lattice design
#**************************
if (method == 'T') {
if (statistic == 'P') {
bound1 = 2^(-0.5) * 3^(-0.25)
bound2 = sqrt(2) * 3^(-0.75)
density = area / ssize
d = sqrt(density)
R = radius / d
if (R <= bound1) pdetect = pi * R^2
if ((R > bound1) & (R <= bound2)) {
Z = acos(bound1 / R)
pdetect = pi*R^2 - 6*Z*R^2 + sqrt(6*sqrt(3)*R^2 - 3) }
if (R > bound2) pdetect = 1
distance = 1.075 * d
}

if (statistic == 'N') {
bound1 = 2^(-0.5) * 3^(-0.25)
bound2 = sqrt(2) * 3^(-0.75)
pdetect.1 = pi * bound1^2
if (pdetect <= pdetect.1) R = sqrt(pdetect / pi)
if ((pdetect > pdetect.1) & (pdetect <= 1)) R = optimize(fT.detect, interval=c(bound1, bound2), pdet = pdetect, maximum=F)$min
if (pdetect == 1) R = bound2
#radius=1.1
d = radius / R
#area = 200
ssize = area / d^2
distance = 1.075 * d
}
if (statistic == 'R') {
bound1 = 2^(-0.5) * 3^(-0.25)
bound2 = sqrt(2) * 3^(-0.75)
pdetect.1 = pi * bound1^2
if (pdetect <= pdetect.1) R = sqrt(pdetect / pi)
if ((pdetect > pdetect.1) & (pdetect <= 1)) R = optimize(fT.detect, interval=c(bound1, bound2), pdet = pdetect, maximum=F)$min
if (pdetect == 1) R = bound2
d = sqrt(area / ssize)
radius = R * d
distance = 1.075 * d
}
}
list(prob = pdetect, ssize = ssize, rad = radius, sep = distance)
}


detect.prop = function(statistic, theta=NA, pdetect=NA, ssize=NA) {
  #*****************************************************************************
  # Written by JB 1/9/16
  # Calculates probability of detection of feature covering an area theta
  # of the survey area. Works for vectors.
  #
  # Statistic = 'P' (probability detection)
  #           = 'N' (sample size)
  #           = 'F' (feature proportion)
  #
  # Output    $prob  = probability feature detected
  #           $ssize = sample size
  #           $prop   = feature proportion
  #****************************************************************************
  #***********************************************
  # Error Checks on parameters
  #***********************************************
  allowedstatistics = c("P", "N", "F")
  if( !(statistic %in% allowedstatistics) )
    stop("Statistic must be of type character and one of P, N or F")
  
  #************************************************************************
  # Check that sufficient parameters have been input
  #************************************************************************
  #if (statistic=="P") {
  #  if (is.na(theta)==T | is.na(ssize)==T) {
  #    stop("theta and ssize arguments needed")
  #  }}
  
  # if (statistic=="N") {
  #  if (is.na(theta)==T | is.na(pdetect)==T) {
  #    stop("theta and pdetect arguments needed")
  #  }}
  
  #if (statistic=="F") {
  #  if (is.na(ssize)==T | is.na(pdetect)==T) {
  #    stop("ssize and pdetect arguments needed")
  #  }}
  
  #************************************************************************
  
  if (statistic == "N") {
    for (j in 1:length(pdetect)) {
      if ( pdetect[j] <=0 | pdetect[j] >=1)
        stop("pdetect must be greater than 0 and less than 1")
    }
    for (j in 1:length(theta)) {
      if ( theta[j] <=0 | theta[j] >=1)
        stop("theta must be greater than 0 and less than 1")
    }
    if (length(pdetect) != length(theta))
      stop("pdetect and theta must be of the same length")
  }
  
  if (statistic == "F") {
    for (j in 1:length(pdetect)) {
      if ( pdetect[j] <=0 | pdetect[j] >=1)
        stop("pdetect must be greater than 0 and less than 1")
    }
    for (j in 1:length(ssize)) {
      if ( !is.wholenumber(ssize[j]) | ssize[j]<0.5 )
        stop("ssize must be a positive integer")
    }
    if (length(pdetect) != length(ssize))
      stop("pdetect and ssize must be of the same length")
  }
  
  if (statistic == "P") {
    for (j in 1:length(ssize)) {
      if ( !is.wholenumber(ssize[j]) | ssize[j]<0.5 )
        stop("ssize must be a positive integer")
    }
    for (j in 1:length(theta)) {
      if ( theta[j] <=0 | theta[j] >=1)
        stop("theta must be greater than 0 and less than 1")
    }
    if (length(ssize) != length(theta))
      stop("ssize and theta must be of the same length")
  }
  
  
  #**************************
  # Do the calculations
  #**************************
  if (statistic == 'P') {
    pdetect = 1 - (1 - theta)^ssize
  }
  if (statistic == 'N') {
    ssize = log(1-pdetect) / log(1-theta)
  }
  if (statistic == 'F') {
    theta = 1 - (1-pdetect)^(1/ssize)
  }
  
  list(prob = pdetect, ssize = ssize, prop = theta)
}


