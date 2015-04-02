\name{permute.groups}
\alias{permute.groups}
\title{
Does randomisation test for the difference in means of two vectors \code{v1} and \code{v2}.
}
\description{
Does randomisation test for the difference in means mu1, mu2
of two vectors \code{v1} and \code{v2}. Can do one or two sided tests.
}
\usage{
permute.groups(v1, v2, alternative, nreps)
}
\arguments{
  \item{v1}{Data vector for variable 1} 
  \item{v2}{Data vector for variable 2}
  \item{alternative}{A character string specifying the alternative
hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} (\code{mu1>mu2}) or \code{"less"}.
You can specify just the initial letter.}
  \item{nreps}{Number of replications used in the randomisation and generation of
the p-value. Default is \code{nreps=999}}
}

\details{Under the null hypothesis that \code{mu1=mu2}, the labelling of the \code{n1+n2} observations is unimportant.
Therefore, we can generate the null distribution for the test statistic \code{m1-m2} or \code{|m1-m2|} depending
on whether a one
or two sided test is required) by randomly permuting the treatment labels nreps times and calculating the test statistic
each time. The p-value is calculated as suggested by Manly (2006).
}

%% There should be a better way to do this (see ?t.test)
\value{The p-value is returned as \code{$p.value}}

\references{
Manly BFJ (2006) Randomization, Bootstrap And Monte Carlo Methods in Biology: 3rd edition. Chapman and Hall.}

\author{
Jon Barry: Jon.Barry@cefas.co.uk
}
\seealso{\code{\link{power.groups}}, \code{\link{permute.BACI}}   }

\examples{
set.seed(5)
v1 = rnorm(27,10,2); v2=rnorm(25,11,2)
permute.groups(v1, v2, alternative="two")
permute.groups(v1, v2, alt="l")
}
\keyword{permutation test}

