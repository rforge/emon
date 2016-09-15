\name{GVFCMOM}

\alias{GVFCMOM}

\title{Function to calculate the method of moments visual fast count estimator}

\description{
The function takes data in the form of counts per segment along a transect and uses the 
raw Generalised Visual Fast Count estimator (as calculated by \code{\link{GVFCMOM}} and its
expectation (as calculated by \code{\link{expected.pois}} for Poisson or \code{\link{expected.nb}} for
negative binomial) to calculate
a method of moments estimator. This effectively,
adjusts the biased raw GVFC estimate. The function allows counts to have either a Poisson or a Negative Binomial
distribution. The method is a generalisation of the methods in Barry and Coggan (2010).
}

\usage{
GVFCMOM(counts, s, d, method, k=1, lowint=0, highint=100)
}

\arguments{
  \item{counts}{Vector of length s that contains a count for each segment}
  \item{s}{Number of segments}
  \item{d}{Number of positive segment counts needed before counting stops}
  \item{method}{Whether Poisson (\code{"pois"}) or Negative Binomial Distribution (\code{"nb"}) is assumed}
  \item{k}{Size parameter of the Negative Binomial distribution}
  \item{lowint}{Minimum value for MOM estimate (default=0)}
  \item{highint}{Maximum value for MOM estimate (default=100)}
}

\value{The method of moments estimate for the transect is returned}

\references{
Barry J, Eggleton J, Ware S and Curtis M (2015) Generalizing Visual Fast Count Estimators for
Underwater Video Surveys. Ecosphere. 
http://www.esajournals.org/doi/full/10.1890/ES15-00093.1
}
\author{
Jon Barry: Jon.Barry@cefas.co.uk
}

\examples{
counts = c(0, 0, 0, 0, 1, 1, 1, 2, 1)
GVFCMOM(counts, s=9, d=2, method='nb', lowint=0, highint=100)
GVFCMOM(counts, s=9, d=1, method='nb', lowint=0, highint=100)
GVFCMOM(counts, s=9, d=1, method='pois', lowint=0, highint=100)
}
\keyword{Visual Fast Count}

\seealso{ \code{\link{GVFC}}, \code{\link{expected.pois}}, \code{\link{expected.nb}} }
