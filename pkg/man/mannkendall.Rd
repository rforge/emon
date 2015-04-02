\name{mannkendall}

\alias{mannkendall}

\title{ Mann-kendall test for trend }

\description{
Calculates the Mann-Kendall statistic for monotonic trend and also the p-value against the null hypothesis of no trend.
Unlike the function \code{\link[Kendall]{MannKendall}}, works for repeat values of time.
}

\usage{ mannkendall(time, Y, nsims.mk=999) }

\arguments{
  \item{time}{Vector of values which define the direction of the trend.}

  \item{Y}{Vector of values for which you want to determine the trend.}

 \item{nsims.mk}{Number of replicate permuations to calculate the p-value. Default=999.}
}

\details{
Error checks of parameters not included so as not to slow down mannkendall. The statistic is calculated by considering
each case \code{j} and considering the subset of observations that have time greater than \code{time[j]}. The Mann Kendall
statistic is the number of observations in this subset for which \code{Y > Y[j]} minus the number for
which \code{Y < Y[j]}.
The statistic is summed over all \code{j}. The p-value is calculated by \code{nreps} random permutations of the Y values.
}

\value{
\item{mann}{Mann-Kendall statistic of observed data, as calculated by \code{\link{mannkendall.stat}}.}

\item{pvalue}{P-value assuming a null hypothesis of no trend and two-way alternative hypothesis.}
}

\references{
Mann, H.B. (1945), Nonparametric tests against trend, Econometrica, 13,
245-259.

Kendall, M.G. 1975. Rank Correlation Methods, 4th edition, Charles Griffin, 
London.
}

\author{
Jon Barry: Jon.Barry@cefas.co.uk
}

\seealso{\code{\link{mannkendall.stat}}, \code{\link{power.trend}}, \code{\link[Kendall]{MannKendall}} }

\examples{
x = rep(1:10,rep(2,10))
y = rnorm(20, 5, 2)
mannkendall(x, y)
}

