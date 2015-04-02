\name{expected.pois}
\alias{expected.pois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expected value of Visual Fast Count Estimator assuming Poisson distribution for counts
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function is used to obtain the method of moments estimator within the function \code{\link{GVFCMOM}}.
Calculates the expected value of the Visual Fast Count method.
The function assumes that the count per segment is Poisson
with mean m/s and that segment counts are independent. The expected value is also a function of the number of
positives d before segment counting is stopped.
}

\usage{
expected.pois(m, s, d) 
}

\arguments{
  \item{m}{The underlying mean of the Poisson process per transect}
  \item{s}{The number of segments per transect}
  \item{d}{The numbner of positive counts before segment counting is stopped}
}
\value{The expected value of the Visual Fast Count estimator}

\author{Jon Barry: Jon.Barry@cefas.co.uk}

\seealso{\code{\link{expected.nb}}, \code{\link{GVFCMOM}}}