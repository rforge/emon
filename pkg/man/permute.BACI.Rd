\name{permute.BACI}
\alias{permute.BACI}
\title{
Does non-parametric randomisation test for the interaction term in a BACI design. 
}
\description{
We have control and treatment data from time 1 in a BACI design, plus control and treatment
data from time 2. The interaction is the amount that the difference in the control and
treatment means is different between times 1 and 2.
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
There are several permutations that can be used to generate the null distribution for
the interaction (see Manly, 2006 and Anderson and Ter Braak, 2003). The method used here is
to do a complete randomisation of the raw data.

The p-value is calculated as suggested by Manly (2006).
}

\value{The p-value is returned as \code{$p.value}}

\references{
Manly BFJ (2006) Randomization, Bootstrap And Monte Carlo Methods in Biology: 3rd edition. Chapman and Hall.

Anderson, M.J. and Ter Braak, C.J.F. (2003). Permutation tests for multi-factorial analysis of
variance. Journal of Computation and Simulation, 73, 85-113.}

\author{
Jon Barry: Jon.Barry@cefas.co.uk
}
\seealso{\code{\link{power.BACI}}, \code{\link{permute.groups}}}

\examples{
}
\keyword{permutation test}

