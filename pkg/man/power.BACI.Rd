\name{power.BACI}

\alias{power.BACI}

\title{
Calculates power for a Before and After Control Impact (BACI) design.
}

\description{
BACI designs are commonly used in environmental monitoring. They are relevant where you
want to measure the effect of an impact (e.g. the effect on benthic ecology of dredging
in an area). Observations for treatment and control areas are measured BEFORE and
after the impact. This function allows you to examine the power of particular BACI designs
to detect differences between the control and the treatment.
}

\usage{ power.BACI(change, change.type="M", nt, nc, parst, parsc,
          distribution, test="P", alpha=0.05, nsims=1000)
}
\arguments{

\item{change}{AFTER treatment mean minus BEFORE treatment mean or percentage change of
   AFTER treatment mean relative to BEFORE treatment mean (depending on value of
   \code{change.type}).
}

\item{change.type}{Whether the parameter \code{change} represents an additive ("A") or
     percentage ("M") change.
}

\item{nt}{Vector of sample sizes for treatment group. Must be of same dimension as \code{nc}.
}


\item{nc}{Vector of sample sizes for control group. Must be of same dimension as \code{nt}.
}


\item{parst}{
    Parameters for the treatment data.
   
    If \code{distribution="Normal"}, \code{parst[1]}
    contains the BEFORE mean and \code{parst[2]} contains the BEFORE standard deviation. If
    \code{distribution="Poisson"}, \code{parst[1]} contains the BEFORE mean. If
   \code{distribution="Lognormal"}, \code{parst[1]} contains the BEFORE mean of the natural
    log data and
    \code{parst[2]} contains the BEFORE standard deviation of the log data. If
    \code{distribution="Negbin"},
    \code{parst[1]} contains the BEFORE mean, \code{parst[2]} contains the BEFORE size, and
    \code{parst[3]} is the AFTER size.
}
\item{parsc}{
    Parameters for the control data.
   
    If \code{distribution="Normal"}, \code{parsc[1]}
    contains the BEFORE mean and \code{parsc[2]} contains the BEFORE standard deviation. If
    \code{distribution="Poisson"}, \code{parsc[1]} contains the BEFORE mean. If
    \code{distribution="Lognormal"}, \code{parsc[1]} contains the BEFORE mean of the natural log
    data and
    \code{parsc[2]} contains the BEFORE standard deviation of the log data. If
    \code{distribution="Negbin"},
    \code{parsc[1]} contains the BEFORE mean, and \code{parsc[2]} contains the BEFORE size.
}
  \item{distribution}{
    The statistical distribution for the two groups. Can be either: \code{"Normal"},
    \code{"Poisson"}, \code{"Lognormal"} or \code{"Negbin"}.
}
  \item{test}{The statistical test used to compare the interaction between the control and treatment
   means at the BEFORE and AFTER sampling occasions. If \code{test="NP"}, then the
   test will be a non-parametric randomisation test, in the spirit of Manly (1997), using the function
   \code{\link{permute.BACI}}.

   If \code{test="P"}, then parametric tests are made to
   compare the treatment (i.e. a factor indicating whether an observation is from the treatmment
   or the control) by time (i.e. a factor indicating whether observations are BEFORE or AFTER)
   interaction. If \code{distribution="Normal"} then this is calculated from the usual Analysis
   of Variance. The same method is used (but on the log data) if \code{distribution="Lognormal"}.

   When \code{distribution="Poisson"} or \code{distribution="Negbin"}, interactions from analysis
   of deviance tables are used to measure the interaction. The p-value is calculated by assuming
   that this interaction deviance has a chi-squared distribution on 1 df. For the Negative
   Binomial distribution, terms in the analysis of deviance table use the same
   value for the size parameter as that estimated from the null model.
}
  \item{alpha}{ If the p-value for the interaction is less than \code{alpha}, a change is deemed
   to have been detected. Used to assess power from the \code{nreps} simulations.
}
  \item{nsims}{ Number of repeat simulations used to calculate the power. Default is 1000.}
}
\details{
BACI designs are relevant where you want to measure the effect of an
impact (e.g. the effect on benthic ecology of dredging in an area).
You take a number of samples both in (treatment) and outside (control) the affected area BEFORE the impact. Those in
the control area should be as similar to the treatment area as possible in terms of benthic ecology. You then sample
the areas again AFTER the impact has taken place. If there is an interaction for the 2x2 crossed design then there is
an effect of the impact. That is, if the control area has changed differently to the treatment area.

This function allows you to examine the power of particular BACI designs. The distribution of the measure being used
can be Normal, Poisson, Negative Binomial or Lognormal.

It is also assumed that the sample sizes before and after the impact are the same, although the sample size in the
treatment area can be different to the control area. Thus, if 10 and 8 samples are taken in the treatment and control
areas before the impact, then 10 and 8 samples are assumed to be taken after the impact.

For the Normal, Poisson and Negative Binomial distributions, the parameter change is simply the percentage or
additive change of the treatment mean from
the BEFORE to the AFTER sampling occasions. For the Lognormal distribution, the percentage change is relative to the
BEFORE treatment mean on the non-log scale. The BEFORE treatment mean is estimated from the mean of the
log data (\code{parst[1]}), the standard deviation of the log data (\code{parst[2]}) and the proposed sampling size
\code{nt}.

The estimator used is the one proposed by Shen (2006), which did better in terms of mean squared error than both the
sample mean on the non-log scale and the maximum likelihood estimators. This is given by
\code{mean.before = exp(parst[1] + (nt-1)*ss / (2*(nt+4)*(nt-1) + 3*ss)), where ss = (nt-1)*parst[2]**2}.

The Negative Binomial distribution option (\code{parst[3]} allows the user to specify the size parameter of the
AFTER treatment distribution. One possibility is to keep the size the same for both the BEFORE and AFTER
distributions. However, because the mean changes and because the variance V = mu+mu^2/size, this means that V
will be different for the BEFORE and AFTER distributions. If you want to keep the variance the same, you can use
the function \code{\link{size2.samevar}}.

Several powers can be calculated per call to this function by specifying more than one values for
the sample sizes \code{nt} and \code{nc}.
}
\value{
  \item{power}{ The estimated power for the design. }

  \item{before.mean}{ The treatment mean used for the before sampling. Only really of interest if
                         \code{method="Lognormal"} as this gives the Shen estimator.}
}

\references{
Shen H, Brown LD and Zhi H (2006) Efficient estimation of log-normal means with application to pharmacokinetic data.
Statistics in Medicine, 25, 3023 to 3038.
}

\author{
Jon Barry (email Jon.Barry@cefas.co.uk)
}

\seealso{
\code{\link{power.trend}}, \code{\link{power.groups}}, \code{\link{size2.samevar}}
}

\examples{

# Data is richness (number of species) and abundance from grab samples from
# the Dogger Bank, UK.

rich = c(15,16,37,12,15,5,13,16,17,34,23,20,22,30,85,55,13,19,30,41,22,8,43,10,38,24,17,
         23,17,17,24,33,30,18,26,18,12,50,19,21,35)
abun = c(50,91,140,21,25,8,28,37,30,90,56,50,40,83,964,180,21,60,81,138,67,17,250,63,152,
         68,42,69,57,67,74,96,75,44,61,49,62,281,55,50,198)

par(mfrow=c(2,2))
hist(rich)
hist(abun)
hist(sqrt(rich))
hist(log(abun))

ssize = seq(10, 50, 10)
parsc.rich = mean(rich); parst.rich = mean(rich)
parsc.abun = rep(0,2); parst.abun = parsc.abun
parst.abun[1] = mean(log(abun)); parst.abun[2] = sd(log(abun))
parsc.abun = parst.abun
power.rich = rep(0, length(ssize))
power.abun = rep(0, length(ssize))

power.rich = power.BACI(change=35, change.type="M", nt=ssize, nc=ssize,
                parst=parst.rich, parsc=parsc.rich,
                distribution="Poisson", test="P", nsims=500)$power

power.abun = power.BACI(change=35, change.type="M", nt=ssize, nc=ssize,
                parst=parst.abun, parsc=parsc.abun, distribution="Lognormal",
                test="P", nsims=500)$power

par(mfrow=c(1,1))
plot(ssize, power.rich, ylim=c(0,1), ylab="Power", xlab="Sample size", type="l")
lines(ssize, power.abun, lty=2, col=2)
legend("bottomright", legend=c("Richness power", "Abundance power"), lty=c(1,2),
       col=c(1,2))
title("BACI power plots")
}