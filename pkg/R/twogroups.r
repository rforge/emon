permute.groups = function(v1, v2, alternative="t", nreps=999) {
#*************************************************************
# Does randomisation test for the difference in means 
# of two vectors v1 and v2.
# alternative can be "two-sided", "greater" or "less"
#*************************************************************
alt = substring(alternative,1,1)

allowedalts = c("l", "t", "g")
if( !(alt %in% allowedalts) )
      stop("Tests must be of type character and one of l, t or g")

if (!is.numeric(v1)) stop("v1 must be a numeric vector")
if (!is.numeric(v2)) stop("v2 must be a numeric vector")

if ( !is.wholenumber(nreps) | nreps<0.5 )
      stop("nreps must be a positive integer")

v1 = v1[!is.na(v1)]
v2 = v2[!is.na(v2)]

obs.diff = mean(v1) - mean(v2)
v12 = c(v1,v2)
l12 = length(v12)
l1 = length(v1) 
sim.diff = rep(0,nreps)
for (j in 1:nreps) {
perm = sample(v12)
sim.diff[j] = mean(perm[1:l1]) - mean(perm[(l1+1):l12])
}
if (alt=="g") bigger = sim.diff[sim.diff >= obs.diff]
if (alt=="l") bigger = sim.diff[sim.diff <= obs.diff]
if (alt=="t") bigger = sim.diff[abs(sim.diff) >= abs(obs.diff)]

pvalue = (length(bigger)+1)/(nreps+1)
list(p.value=pvalue)
}


power.groups = function(change, change.type="M", n1, n2, pars1, pars2,
           distribution, test, alternative="two", alpha=0.05, nsims=1000,
           nreps=999)  {
#********************************************************************
# distribution = "Negbin" needs library MASS
#********************************************************************

alloweddistributions = c("Normal", "Poisson", "Negbin", "Lognormal")
if( !(distribution %in% alloweddistributions) )
      stop("Distribution must be of type character and be one of Normal, Poisson, Negbin or Lognormal")

allowedtests = c("P", "NP")
if( !(test %in% allowedtests) )
      stop("Tests must be of type character and be one of P or NP")

allowedchange.types = c("A", "M")
if( !(change.type %in% allowedchange.types) )
      stop("Change.type must be of type character and be one of A or M")

length.n1 = length(n1); length.n2 = length(n2)

if ( length.n1 != length.n2 ) stop("n1 and n2 must be of the same length")

for (j in 1: length.n1) {
if ( !is.wholenumber(n1[j]) | n1[j]<0.5 )
      stop("n1 must contain positive integers")

if ( !is.wholenumber(n2[j]) | n2[j]<0.5 )
      stop("n2 must contain positive integers")
}

if (distribution=="Normal") {
    if (length(pars1)!=2) stop("pars1 must be of length 2")
    if (length(pars2)!=1) stop("pars2 must be of length 1")
    if (pars1[2]<=0) stop("pars1[2] must be > 0")
    if (pars2[1]<=0) stop("pars2[1] must be > 0")
    varequal = F
    if (pars1[2] == pars2[1]) varequal = T
}

if (distribution=="Poisson") {
    if (length(pars1)!=1) stop("pars1 must be of length 1")
    if (pars1[1]<0) stop("pars1[1] must be > 0")
}

if (distribution=="Lognormal") {
    if (length(pars1)!=2) stop("pars1 must be of length 2 for Lognormal")
    if (length(pars2)!=1) stop("pars2 must be of length 1 for Lognormal")
    if (pars1[2]<=0) stop("Group 1 standard deviation must be positive")
    if (pars2[1]<=0) stop("Group 2 standard deviation must be positive")
    varequal = F
    if (pars1[2] == pars2[1]) varequal = T
}

if (distribution=="Negbin") {
    if (length(pars1)!=2) stop("pars1 must be of length 2 for Negative Binomial")
    if (length(pars2)!=1) stop("pars2 must be of length 1 for Negative Binomial")
    if (pars1[2]<=0) stop("Group 1 standard deviation must be positive")
    if (pars2[1]<=0) stop("Group 2 standard deviation must be positive")
}

alt = substring(alternative, 1, 1)
allowedalts = c("l", "t", "g")
if( !(alt %in% allowedalts) )
      stop("Tests must be of type character and one of l, t or g")

if (alpha <= 0 | alpha >= 1) stop("Type 1 error must be between 0 and 1")

if ( !is.wholenumber(nreps) | nreps<0.5 )
      stop("nreps must be a positive integer")

power = vector("numeric", length.n1)

if (distribution=="Normal") {
m1 = pars1[1]
if (change.type=="M") m2 = m1 + change*m1/100
if (change.type=="A") m2 = m1 + change
s1 = pars1[2]
s2 = pars2[1]
for (k in 1:length.n1) {
count = 0
for (j in 1:nsims) {
y1 = rnorm(n1[k], m1, s1)
y2 = rnorm(n2[k], m2, s2)
if (test=='P') p = t.test(y1, y2, var.equal=varequal, alternative=alt)$p.val
if (test=='NP') p = permute.groups(y1, y2, alternative=alt, nreps=nreps)$p.val
if (p < alpha) count = count + 1 }
power[k] = count / nsims
}}

if (distribution=="Poisson") {
m1 = pars1[1]
if (change.type=="M") m2 = m1 + change*m1/100
if (change.type=="A") m2 = m1 + change
for (k in 1:length.n1) {
x = factor(c(rep(1, n1[k]), rep(2, n2[k])))
count = 0
for (j in 1:nsims) {
y1 = rpois(n1[k], m1)
y2 = rpois(n2[k], m2)
if (test=='NP') p = permute.groups(y1, y2, alternative=alt, nreps=nreps)$p.val
if (test=='P') {
y = c(y1, y2)
dev0 = glm(y ~ 1, family=poisson)$dev
dev1 = glm(y ~ x, family=poisson)$dev
p = (1 - pchisq(dev0-dev1, 1))
}
if (p < alpha) count = count + 1 }
power[k] = count / nsims
}}

if (distribution=="Negbin") {
m1 = pars1[1]
if (change.type=="A") m2 = m1 + change
if (change.type=="M") m2 = m1 + change*m1/100
size1 = pars1[2]
size2 = pars2[1]
for (k in 1:length.n1) {
x = factor(c(rep(1, n1[k]), rep(2, n2[k])))
count = 0
for (j in 1:nsims) {
y1 = rnbinom(n=n1[k], size=size1, mu=m1)
y2 = rnbinom(n=n2[k], size=size2, mu=m2)
if (test=='NP') p = permute.groups(y1, y2, alternative=alt, nreps=nreps)$p.val
if (test=='P') {
y = c(y1, y2)
temp = glm.nb(y ~ x)
p = anova(temp)[2,5]
}
if (p < alpha) count = count + 1 }
power[k] = count / nsims
}}

if (distribution=="Lognormal") {
# mu and sigma are on log scale
mu1 = pars1[1]
sigma1 = pars1[2]
sigma2 = pars2[1]
mean1 = exp(mu1 + (sigma1^2) / 2)
if (change.type=="M") mean2 = mean1 + change*mean1/100
if (change.type=="A") mean2 = mean1 + change
mu2 = log(mean2) - (sigma2^2)/2
for (k in 1:length.n1) {
count = 0
for (j in 1:nsims) {
y1 = rlnorm(n1, mu1, sigma1)
y2 = rlnorm(n2, mu2, sigma2)
if (test=='NP') p = permute.groups(y1, y2, alternative=alt, nreps=nreps)$p.val
if (test=='P') p = t.test(log(y1), log(y2), var.equal=varequal, alternative=alt)$p.val
if (p < alpha) count = count + 1 }
power[k] = count / nsims
}}
power
}


size2.samevar = function(mu1, mu2, s1) {
#********************************************************************
# Calculates the size for group 2 such that the Negative
# Binomial variance is the same in the two groups. This value can
# then be fed into power.groups.
#********************************************************************
s2 = mu2^2 * s1 / (s1*(mu1 - mu2) + mu1^2)
if (s2 <=0 | s2==Inf) stop("No positive, finite solution for size2")
s2
}

