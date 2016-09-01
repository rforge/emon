is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

power.BACI = function(change, change.type="M", nt, nc, parst, parsc,
          distribution, test="P", alpha=0.05, nsims=1000)  {
#********************************************************************
# Distribution can be Normal, Poisson Negative Binomial or Lognormal
# nc = n for control
# nt = n for treatment
# Assume same sample sizes before and after
# For Lognormal, the first two entries of parst() and parsc() are on the
# log scale. The third, for parst() is the mean on the original scale.
#
# The % reduction is on the original scale?
#
# Have I got the variance right when it is changed Think OK as I am
# assuming no change on the log scale.
#********************************************************************
alloweddistributions = c("Normal", "Poisson", "Negbin", "Lognormal")
if( !(distribution %in% alloweddistributions) )
      stop("Distribution must be one of Normal, Poisson, Negbin or
            Lognormal")

allowedchange.types = c("A", "M")
if( !(change.type %in% allowedchange.types) )
      stop("Change.type must be of type character and be one of A or M")

allowedtests = c("P", "NP")
if( !(test %in% allowedtests) )
      stop("Tests must be of type character and be one of P or NP")

if (distribution=="Normal") {
    if (length(parst)!=2) stop("parst must be of length 2")
    if (length(parsc)!=2) stop("parsc must be of length 2")
    if (parst[2]<=0) stop("Treatment standard deviation must be positive")
    if (parsc[2]<=0) stop("Control standard deviation must be positive")
}

if (distribution=="Poisson") {
    if (length(parst)!=1) stop("parst must be of length 1")
    if (length(parsc)!=1) stop("parsc must be of length 1")
    if (parst[1]<=0) stop("Treatment mean must be positive")
    if (parsc[1]<=0) stop("Control mean must be positive")
}

if (distribution=="Negbin") {
    if (length(parst)!=3) stop("parst must be of length 3")
    if (length(parsc)!=2) stop("parsc must be of length 2")
    if (parst[1]<=0) stop("Treatment mean must be positive")
    if (parsc[1]<=0) stop("Control mean must be positive")
    if (parst[2]<=0) stop("Treatment size must be positive")
    if (parsc[2]<=0) stop("Control size must be positive")
    if (parst[3]<=0) stop("Treatment size (2) must be positive")
}

if (distribution=="Lognormal") {
    if (length(parst)!=2) stop("parst must be of length 2")
    if (length(parsc)!=2) stop("parsc must be of length 2")
    if (parst[2]<=0) stop("Treatment standard deviation must be positive")
    if (parsc[2]<=0) stop("Control standard deviation must be positive")
}

if (alpha<0 | alpha>1) stop("alpha must lie between 0 and 1")

if ( !is.wholenumber(nsims) | nsims<0.5 )
      stop("nsims must be a positive integer")

length.nt = length(nt); length.nc = length(nc)

if ( length.nt != length.nc ) stop("nt and nc must be of the same length")

for (j in 1: length.nt) {
if ( !is.wholenumber(nt[j]) | nt[j]<1.5 )
      stop("nt must contain positive integers")

if ( !is.wholenumber(nc[j]) | nc[j]<1.5 )
      stop("nc must contain positive integers greater than 1")
}

#### END OF ERROR-CHECKING STUFF ####

power = vector("numeric", length.nt)

if (distribution=="Normal") {
before.treat.mean = parst[1]
mu.treat = parst[1]; sd.treat = parst[2]
mu.con = parsc[1]; sd.con = parsc[2]
if (change.type=="M") mu.treat.after = mu.treat + change * mu.treat / 100
if (change.type=="A") mu.treat.after = mu.treat + change

for (k in 1:length.nt) {
count = 0
group = factor(c(rep(1,nt[k]), rep(2,nc[k]), rep(1,nt[k]), rep(2,nc[k])))
time = factor(c(rep(1,nt[k]), rep(1,nc[k]), rep(2,nt[k]), rep(2,nc[k])))

for (j in 1:nsims) {
con.before = rnorm(nc[k], mu.con, sd.con)
treat.before = rnorm(nt[k], mu.treat, sd.treat)
con.after = rnorm(nc[k], mu.con, sd.con)
treat.after = rnorm(nt[k], mu.treat.after, sd.treat)
if (test=="P") {
y = c(treat.before, con.before, treat.after, con.after)
arse = aov(y ~ group*time)
p = anova(arse)[3,5]
}
if (test=="NP") p = permute.BACI(treat.before, con.before, treat.after,
                    con.after)
if (p < alpha) count = count + 1
}
power[k] = count / nsims
}}

if (distribution=="Poisson") {
before.treat.mean = parst[1]
mu.con = parsc[1]; mu.treat = parst[1]
if (change.type=="M") mu.treat.after = mu.treat + change * mu.treat / 100
if (change.type=="A") mu.treat.after = mu.treat + change

for (k in 1:length.nt) {
count = 0
group = factor(c(rep(1,nt[k]), rep(2,nc[k]), rep(1,nt[k]), rep(2,nc[k])))
time = factor(c(rep(1,nt[k]), rep(1,nc[k]), rep(2,nt[k]), rep(2,nc[k])))

for (j in 1:nsims) {
con.before = rpois(nc[k], mu.con)
treat.before = rpois(nt[k], mu.treat)
con.after = rpois(nc[k], mu.con)
treat.after = rpois(nt[k], mu.treat.after)
if (test=="P") {
y = c(treat.before, con.before, treat.after, con.after)
temp = glm(y ~ group*time, family=poisson)
p = 1 - pchisq(anova(temp)[4,2],1)
}
if (test=="NP") p = permute.BACI(treat.before, con.before, treat.after,
                    con.after)
if (p < alpha) count = count + 1
}
power[k] = count / nsims
}}

if (distribution=="Negbin") {
mu.con = parsc[1]
mu.treat = parst[1]
size.con = parsc[2]
size.treat = parst[2]
size.treat.after = parst[3]
if (change.type=="M") mu.treat.after = mu.treat + change * mu.treat / 100
if (change.type=="A") mu.treat.after = mu.treat + change

for (k in 1:length.nt) {
count = 0
group = factor(c(rep(1,nt[k]), rep(2,nc[k]), rep(1,nt[k]), rep(2,nc[k])))
time = factor(c(rep(1,nt[k]), rep(1,nc[k]), rep(2,nt[k]), rep(2,nc[k])))

for (j in 1:nsims) {
con.before = rnbinom(n=nc[k], size=size.con, mu=mu.con)
treat.before = rnbinom(n=nt[k], size=size.treat, mu=mu.treat)
con.after = rnbinom(n=nc[k], size=size.con, mu=mu.con)
treat.after = rnbinom(n=nt[k], size=size.treat.after, mu=mu.treat.after)

if (test=="P") {
y = c(treat.before, con.before, treat.after, con.after)
temp = glm.nb(y ~ group*time)
p = anova(temp)[4,5]
}
if (test=="NP") p = permute.BACI(treat.before, con.before, treat.after,
                    con.after)
if (p < alpha) count = count + 1
}
power[k] = count / nsims
}}

if (distribution=="Lognormal") {
mu.treat = parst[1]; sigma.treat = parst[2]
mu.con = parsc[1]; sigma.con = parsc[2]
mean.treat = exp(mu.treat + (sigma.treat^2) / 2)

if (change.type=="M") mean.treat.after = mean.treat + change*mean.treat / 100
if (change.type=="A") mean.treat.after = mean.treat + change

mu.treat.after = log(mean.treat.after) - (sigma.treat^2)/2

for (k in 1:length.nt) {
count = 0
group = factor(c(rep(1,nt[k]), rep(2,nc[k]), rep(1,nt[k]), rep(2,nc[k])))
time = factor(c(rep(1,nt[k]), rep(1,nc[k]), rep(2,nt[k]), rep(2,nc[k])))

for (j in 1:nsims) {
con.before = rlnorm(nc[k], mu.con, sigma.con)
treat.before = rlnorm(nt[k], mu.treat, sigma.treat)
con.after = rlnorm(nc[k], mu.con, sigma.con)
treat.after = rlnorm(nt[k], mu.treat.after, sigma.treat)
if (test=="P") {
y = log(c(treat.before, con.before, treat.after, con.after))
temp = aov(y ~ group*time)
p = anova(temp)[3,5]
}
if (test=="NP") p = permute.BACI(treat.before, con.before, treat.after,
                    con.after)
if (p < alpha) count = count + 1
}
power[k] = count / nsims
}}
list(power=power)
}




permute.BACI = function(t1, c1, t2, c2, nreps=999) {
#****************************************************************
# Does randomisation test for the interaction in a BACI design.
#****************************************************************
if (!is.numeric(t1)) stop("t1 must be a numeric vector")
if (!is.numeric(c1)) stop("c1 must be a numeric vector")
if (!is.numeric(t2)) stop("t2 must be a numeric vector")
if (!is.numeric(c2)) stop("c2 must be a numeric vector")

if ( !is.wholenumber(nreps) | nreps<0.5 )
      stop("nreps must be a positive integer")

t1 = t1[!is.na(t1)]
t2 = t2[!is.na(t2)]
c1 = c1[!is.na(c1)]
c2 = c2[!is.na(c2)]

obs.diff = (mean(t1) - mean(c1)) - (mean(t2) - mean(c2))
treat = c(t1,t2)
ltreat = length(treat)
control = c(c1,c2)
lcontrol = length(control)
lt1 = length(t1)
lc1 = length(c1)

sim.diff = rep(0,nreps)
for (j in 1:nreps) {
permt = sample(treat)
permc = sample(control)
sim.diff[j] = (mean(permt[1:lt1]) - mean(permc[1:lc1])) - (mean(permt[(lt1+1):ltreat]) - mean(permc[(lc1+1):lcontrol]))
}
bigger = sim.diff[abs(sim.diff) >= abs(obs.diff)]

pvalue = (length(bigger)+1)/(nreps+1)
list(p.value=pvalue)
}



