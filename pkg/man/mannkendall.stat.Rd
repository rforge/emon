\name{mannkendall.stat}
\alias{mannkendall.stat}
\title{
Mann-Kendall statistic.
}
\description{
Calculates the Mann-Kendall statistic for monotonic trend. Unlike the function \code{\link[Kendall]{MannKendall}}, works
for repeat values of time. Used in function \code{\link{mannkendall}}, which also calculates the p-value by simulation.
}

\usage{ mannkendall.stat(time, Y) }

\arguments{
  \item{time}{ Vector of values which define the direction of the trend. }
  \item{Y}{Vector of values for which you want to determine the trend.}
}

\details{
The statistic is calculated by considering
each case \code{j} and considering the subset of observations that have time greater than \code{time[j]}. The Mann Kendall
statistic is the number of observations in this subset for which \code{Y > Y[j]} minus the number for
which \code{Y < Y[j]}.
The statistic is summed over all \code{j}. The p-value is calculated by \code{nreps} random permutations of the Y values.
}

\value{Mann-Kendall statistic}

\references{
Mann, H.B. (1945), Nonparametric tests against trend, Econometrica, 13,
245-259.
Kendall, M.G. 1975. Rank Correlation Methods, 4th edition, Charles Griffin, 
London.
}
\author{
Jon Barry: Jon.Barry@cefas.co.uk
}
\seealso{\code{\link{mannkendall}}, \code{\link{power.trend}}, \code{\link[Kendall]{MannKendall}}}
