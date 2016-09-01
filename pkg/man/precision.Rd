\name{precision}
\alias{precision}
\title{
Sample size for given precision or precision for given sample size
}
\description{
Precision is measured by the width of a 100(1-alpha)% confidence interval for a single population mean.
The function generates the sample size needed to achieve this or the precision achieved for
a specified sample size.
}
\usage{
precision(d, n, pars, method="sample size", alpha=0.05, minint=1, maxint=500) 
}

\arguments{
\item{d}{
The Confidence Interval width required (for use with \code{method="sample size"}). This can
be a vector.
}
\item{n}{
Sample size (for use with \code{method="width"}). This can be a vector.
}
\item{pars}{
Standard deviation of the variable
}
\item{method}{
 Whether sample size is required (\code{"sample size"}) or precision (\code{"width"}).
}
\item{alpha}{
Defines the (1-alpha/2) percentage point of the t-dristribution used in the confidence interval.
}
\item{minint}{
Lower bound to be used in the search interval for the sample size.
}
\item{maxint}{
Upper bound to be used in the search interval for the sample size.
}
}


\details{
The width of a Confidence Interval for the mean is given by the standard formula
\code{d = 2 * sigma * t(1-alpha/2, n-1) / sqrt(n)}, where sigma is the standard deviation and
n is the sample
size. t(.) is the relevant quantile of the t distribution function.If sample size is required then
we can turn this equation round to get \code{n = [2 * sigma * t(1-alpha/2, n-1)/d]^2}. To solve this
equation for the sample size \code{n}, \code{precision} uses the function \code{optimize}.
}

\value{
\item{n}{
Sample sizes.
}
\item{d}{
Confidence interval widths.
}
}

\author{
Jon Barry: Jon.Barry@cefas.co.uk
}

\examples{

precision(d=c(1,1.2,1.5), pars=1, method="sample size", alpha=0.05)

precision(d=c(4), pars=1, method="sample size", alpha=0.05)

precision(n=c(20,25), pars=1, method="width", alpha=0.05)

}
