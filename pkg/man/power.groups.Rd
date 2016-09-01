\name{power.groups}
\alias{power.groups}
\title{
Power for comparing mean of two groups
}
\description{
Calculates the power by simulation for comparing the mean of two groups of independent observations.
}
\usage{
power.groups(change, change.type="M", n1, n2, pars1, pars2,
           distribution, test, alternative="two", alpha=0.05, nsims=1000, nreps=999)
}

\arguments{
\item{change}{
Mean of second group minus mean of first group (i.e. \code{mu2-mu1}) or percentage change
 in \code{mu1} to create \code{mu2} (depending on value of \code{change.type}).
}
\item{change.type}{
Whether the parameter \code{change} represents an additive ("A") or percentage ("M") change.
}
\item{n1}{
Vector of sample sizes for group 1. Must be of same dimension as \code{n2}.
}
\item{n2}{
Vector of sample sizes for group 2. Must be of same dimension as \code{n1}.
}

\item{pars1}{
    Parameters for the treatment data. If \code{distribution="Normal"}, \code{pars1[1]}
    contains the mean for group 1 and \code{pars1[2]} contains the standard deviation. If
    \code{distribution="Poisson"}, \code{pars1[1]} contains the mean for group 1. If
    \code{distribution="Lognormal"}, \code{pars1[1]} contains the group 1 mean of the natural log
    data and
    \code{pars1[2]} contains the standard deviation of the log data. If \code{distribution="Negbin"},
    \code{pars1[1]} contains the mean and \code{pars1[2]} contains the group 1 size parameter.
}
\item{pars2}{
    \code{pars2[1]} is the standard deviation for group 2 if \code{distribution="Normal"}.
    If \code{distribution="Lognormal"}, \code{pars2[1]} is the standard deviation of the log data for
    group 2. For \code{distribution="Negbin"}, \code{pars2[1]} gives the group 2 size.
}
\item{distribution}{
    The statistical distribution for the two groups. Can be either: \code{"Normal"},
    \code{"Poisson"}, \code{"Lognormal"} or \code{"Negbin"}.
}
\item{test}{
The statistical test used to compare the group means. If \code{test="NP"} then the
test will be a
non-parametric randomisation test, in the spirit of Manly (1997), using the function
\code{\link{permute.groups}}. If \code{test="P"}, then parametric tests are made to
compare the group means. If \code{distribution="Normal"}, a two sample t-test is
carried out. If the standard deviations (defined by \code{pars1[2]} and \code{pars2[1]}
are equal, then the t-test calculates the usual pooled standard deviation. However,
if the standard deviations are not equal then the default method used by
\code{\link{t.test}} is adopted.

When \code{distribution="Lognormal"}, natural logs are taken of the simulated data and
a t-test used in a similar way as to when \code{distribution="Lognormal"}.
When \code{distribution="Poisson"}, the difference in deviances between the null model
and the model with group membership fitted as factor is compared againt a chi-squared
distributiuon on 1 degree of freedom. A similar (but not quite) method is used for when
\code{distribution="Negbin"}. The Generalised Linear Model function
\code{\link[MASS]{glm.nb}} for the Negative Binomial distribution is used. The p-value
for comparing the two groups is taken from the analysis of deviance table after the
model with group membership is fitted as a factor. This p-value, however, uses the same
value for the size parameter, as estimated from the null model, for group member deviance.
This seems to be the correct thing to do as estimating separate size parameters for the two
models mucks up the nesting of the models.
}
  \item{alternative}{
A character string specifying the alternative hypothesis, must be one of \code{"two.sided"}
(default), \code{"greater"} or \code{"less"}. You can specify just the initial letter. As
an example, \code{"less"} means that the alternative is that the mean of the first group
is less than the mean of the second.
}
  \item{alpha}{
The type 1 error for assessing statistical significance (default is 0.05) in the power simulations.
}
  \item{nsims}{
Number of repeat simulations to estimate power (default is 1000).
}
  \item{nreps}{
Number of repeat permutations for randomisation test (default is 999).
}
}
\details{
The Negative Binomial distribution option allows the user to specify the size parameter for both
groups 1 and 2. One possibility is to keep the size the same for both groupss. However, because the
mean is different between the groups and because the variance V = mu+mu^2/size, this means that V
will be different for the group 1 and group 2 distributions. If you want to keep the variance the
same, you can use the function \code{\link{size2.samevar}}.

Several powers can be calculated per call to this function by specifying more than one values for
the sample sizes \code{n1} and \code{n2}.
}
\value{
The power is returned. This is the proportion of the \code{nreps} simulations that returned
a p-value less than the type1 error.
}
\references{
Manly BFJ (1997) Randomization, bootstrap and monte carlo methods in biology:
2nd edition. Chapman and Hall, London, 399 pages.
}
\author{
Jon Barry: Jon.Barry@cefas.co.uk
}

\seealso{
\code{\link{permute.groups}}, \code{\link[MASS]{glm.nb}}, \code{\link{size2.samevar}}
}
\examples{
library(MASS)

power.groups(change=2.5, change.type="A", n1=20, n2=20, pars1=c(10,2),
       pars2=2, test='P', distribution="Normal")

power.groups(change=2.5, change.type="A", n1=seq(5,25,5), n2=seq(5,25,5), pars1=c(10,2),
       pars2=2, test='P', distribution="Normal")

power.groups(change=25, change.type="M", n1=20, n2=20, pars1=10,
       test='P', distribution="Poisson")

power.groups(change=4, change.type="A", n1=20, n2=20, pars1=c(1,2),
       pars2=2, test='P', distribution="Lognormal")


# Keeping size constant
power.groups(change=100, change.type="M", n1=20, n2=20, pars1=c(5,2),
       pars2=2, test='P', distribution="Negbin")

# Keeping variance constant
s2 = size2.samevar(mu1=5, mu2=10, s1=2)   # 13.333
power.groups(change=100, change.type="M", n1=20, n2=20, pars1=c(5,2),
       pars2=s2, test='P', distribution="Negbin")

}
