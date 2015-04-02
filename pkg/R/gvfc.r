is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


GVFC = function(counts, s, d) {
#**************************************************************************
# Generalised raw visual fast count estimator
# counts   : vector of length s that contains counts for each segment
# s      : number of segments
# d      : order of the visual count
#**************************************************************************
# COUNTS can be any length up to s
# If COUNTS has more than d positives, then ignore the excess
#**************************************************************************
# ERROR CHECKS
#****************************************************
if (s < 0.5 | !is.wholenumber(s)) stop("s must be a positive integer")
if (d < 0.5 | !is.wholenumber(d)) stop("d must be a positive integer")

len.counts = length(counts)
if (len.counts > s) stop("counts is of higher dimension than number of segments")

npos = sum(counts > 0.5)
cum.pos = cumsum(counts > 0.5)

if (npos < d & len.counts < s) stop("Not enough entries in counts")

if (npos == 0) vfc = 0

if (npos != 0 & npos < d) {
vfc = sum(counts) / len.counts
}
if (npos >= d) {
seg.nums = 1:len.counts
seg = seg.nums[cum.pos==d][1]
vfc = cumsum(counts)[seg] / seg
}
lambda = vfc * s
lambda
}


GVFCMOM = function(counts, s, d, method, k=1, lowint=0, highint=100) {
#*************************************************************************
# Calculates the generalised method of moments estimator.
# Updated September 16th 2014
#*************************************************************************
# ERROR CHECKS
#**********************************************
if (s < 0.5 | !is.wholenumber(s)) stop("s must be a positive integer")
if (d < 0.5 | !is.wholenumber(d)) stop("d must be a positive integer")

len.counts = length(counts)
if (len.counts > s)
        stop("counts is of higher dimension than number of segments s")

npos = sum(counts > 0.5)
cum.pos = cumsum(counts > 0.5)

if (npos < d & len.counts < s) stop("Not enough entries in counts")

allowedmethods = c("nb", "pois")
if( !(method %in% allowedmethods) )
            stop("method must be one of nb or pois")

if (method=="nb") {
if (k<0) stop("k must be negative for method=nb")
}

if (lowint < 0) stop("lowint must not be negative")
if (highint < 0) stop("highint must not be negative")

if (lowint >= highint) stop("lowint must be less than highint")
 
#***********************************************

est1 = GVFC(counts, s, d)

if (est1 > 0)
{if (method=='pois') est2 = optimize(mom.min.pois, interval=c(lowint,highint), smin=s, dmin=d, est1=est1)$min
if (method=='nb') est2 = optimize(mom.min.nb, interval=c(lowint,highint), smin=s, dmin=d, est1=est1, kmin=k)$min
if (est2<0) est2=0
}
if (est1==0) est2=0
est2
}



mom.min.pois = function(x, smin, dmin, est1) {
exp = abs(expected.pois(x,smin,dmin) - est1)
exp
}


mom.min.nb = function(x, smin, dmin, est1, kmin) {
exp = abs(expected.nb(kmin,x,smin,dmin) - est1)
exp
}


expected.pois = function(m,s,d) {
#************************************************************************
# Calculates the expected value of the Visual Fast Count method as a function
# of the true mean m. This assumes a Poisson process for individuals.
# This bias uses the general formula and is a function of the number of
# segments s and the number of positives d
#************************************************************************
p0 = exp(-m/s)
# < d positives found
jvec = 0:(d-1)
sum1 = sum(p0^(s-jvec)*(1-p0)^(jvec-1)*choose(s,jvec)*jvec)
exp1 = sum1*m/s

sum2 = 0
for (j in d:s) {
sum2 = sum2 + p0^(j-d)*choose(j-1,d-1)/j
}
exp2 = sum2*d*m*((1-p0)^(d-1))
#
expected = exp1 + exp2
expected
}


expected.nb = function(k,m,s,d) {
#************************************************************************
# Calculates the expected value of the Visual Fast Count method as a function
# of the true mean m, the number of segments s and the number of positives d
# This assumes a Negative Binomial distributions for counts.
#************************************************************************
p0 = (1 + (m/(s*k)))^(-k)
# < d positives found
jvec = 0:(d-1)
sum1 = sum(p0^(s-jvec)*(1-p0)^(jvec-1)*choose(s,jvec)*jvec)
exp1 = sum1*m/s

sum2 = 0
for (j in d:s) {
sum2 = sum2 + p0^(j-d)*choose(j-1,d-1)/j
}
exp2 = sum2*d*m*((1-p0)^(d-1))
#
expected = exp1 + exp2
expected
}


