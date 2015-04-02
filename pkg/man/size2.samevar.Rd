\name{size2.samevar}
\alias{size2.samevar}
\title{
Calculates negative binomial size to preserve constant variance.
}
\description{
Calculates the Negative Binomial size parameter \code{s2} such that the variance of the distribution
with mean \code{mu2} and size \code{s2} is the same as the Negative Binomial distribution with mean
\code{mu1} and size \code{s1}. This can be useful when computing power for a Negative Binomial
distribution in the packages \code{\link{power.groups}} and \code{\link{power.BACI}}.
}
\usage{
size2.samevar(mu1, mu2, s1)
}
\arguments{
  \item{mu1}{Negative Binomial mean for group 1}
  \item{mu2}{Negative Binomial mean for group 2}
  \item{s1}{Negative Binomial size for group 1}
}
\value{
The size for group 1.
}
\author{
Jon Barry: Jon.Barry@cefas.co.uk
}
\seealso{
\code{\link{power.groups}}, \code{\link{power.BACI}}
}
\examples{
mu1=5; mu2=10; s1=3
s2 = size2.samevar(mu1, mu2, s1)
s2
# Check variances are the same
v1 = mu1 + mu1^2 / s1
v2 = mu2 + mu2^2 / s2
v1; v2
}

