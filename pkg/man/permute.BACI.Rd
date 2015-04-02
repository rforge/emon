\name{permute.BACI}
\alias{permute.BACI}
\title{
Does non-parametric randomisation test for the interaction term in a BACI design. 
}
\description{
We have control and treatment data from time 1 in a BACI design, plus control and treatment
data from time 2. The interaction  the amount that the difference in the control and
treatment meansis different between times 1 and 2.
}
\usage{
permute.BACI(t1, c1, t2, c2, nreps=999)
}
\arguments{
  \item{t1}{Data vector for the treatment at time 1} 
  \item{c1}{Data vector for the control at time 1}
  \item{t2}{Data vector for the treatment at time 2} 
  \item{c2}{Data vector for the control at time 2} 
  \item{nreps}{Number of replications used in the randomisation and generation of
              the p-value. Default is \code{nreps=999}}
}

\details{
The test statistic used to define the interaction is T=[mean(t1)-mean(c1)]-[mean(t2)-mean(c2)].
We assume a null hypothesis of a zero interaction. For a zero interaction to occur, we need
mean(t1)-mean(c1) = mean(t2)-mean(c2). Thus, the appropriate permutation is to permute the time
labels for the control observations and, similarly, to permute the time labels for the treatment
observations. This creates a null distribution for the interaction (using sample means to calculate
T).
 
The p-value is calculated as suggested by Manly (2006).
}

\value{The p-value is returned as \code{$p.value}}

\references{
Manly BFJ (2006) Randomization, Bootstrap And Monte Carlo Methods in Biology: 3rd edition. Chapman and Hall.}

\author{
Jon Barry: Jon.Barry@cefas.co.uk
}
\seealso{\code{\link{power.BACI}}, \code{\link{permute.groups}}}

\examples{
set.seed(5)
v1 = rnorm(27,10,2); v2=rnorm(25,11,2)
}
\keyword{permutation test}

