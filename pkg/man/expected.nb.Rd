\name{expected.nb}
\alias{expected.nb}

\title{
Expected value of Visual Fast Count Estimator assuming Negative Binomial distribution for counts
}
\description{
The function is used to obtain the method of moments estimator within the function \code{\link{GVFCMOM}}.

Calculates the expected value of the Visual Fast Count method.
The function assumes that the count per segment is Negative Binomial with mean \code{m/s} and size \code{k},
and that segment counts are independent. The expected value is also a function of the number of positives
\code{d} before segment counting is stopped.
}

\usage{
expected.nb(k, m, s, d)
}

\arguments{
  \item{k}{Size parameter of the Negative Binomial distribution}
  \item{m}{The mean of the Negative Binomial distribution per transect}
  \item{s}{The number of segments per transect}
  \item{d}{The number of positive counts before segment counting is stopped}
}
\value{The expected value of the Visual Fast Count estimator}

\author{Jon Barry: Jon.Barry@cefas.co.uk}

\seealso{\code{\link{expected.pois}}, \code{\link{GVFCMOM}}}

