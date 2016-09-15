\name{svariog}
\alias{svariog}
\title{
Calculates empirical semi-variogram.
}
\description{
Calculates empirical semi-variogram cloud plus classical, robust and median estimators from bins.
}

\usage{
           svariog(x, y, z, u) 
}

\arguments{
  \item{x}{Location vector 1 (e.g. longitude).} 
  \item{y}{Location vector 2 (e.g. latitude).}
  \item{z}{Response vector observed at the locations.}
  \item{u}{\code{(b+1)} cut points used to define the \code{b} bins for distances. The cut points define the
           boundaries
           for each bin.  Distances on the boundary of bins go into the lower bin.}
}

\details{Generates the \code{n(n-1)/2} distances between each of the n points together with the semi-variogram
         cloud of
         the \code{n(n-1)/2} differences \code{(zi-zj)^2 / 2} between pairs of observations \code{(i,j)}.
         This cloud
         is smoothed by taking one of three sorts of averages
         within each bin - bin end points are defined by the vector \code{u}. These averages are the
         classical (the
         bin mean) estimator, a function of the bin median and a robust estimator. Both the median and the
         robust estimators are based on absolute differences between \code{z} pairs. These methods are
         defined in Cressie (1993).}

\value{
\item{classical}{Classical semi-variogram estimator.}

\item{med}{Median semi-variogram estimator.}

\item{robust}{Robust semi-variogram estimator.}

\item{freq}{Frequencies of distances within each bin.}

\item{mid}{Mid points of each bin.}

\item{zcloud}{Unsmoothed semi-variogram cloud.}

\item{dcloud}{Distances between pairs of points for the variogram cloud.} 
}

\references{
Cressie, NAC (1993) Statistics for Spatial Data, Revised Edition. Wiley, New York.}

\author{
Jon Barry: Jon.Barry@cefas.co.uk
}
\seealso{\code{\link[geoR]{variog}}}

\examples{
# Example based on the number of benthic species found from samples of Hamon Grabs from 50
# locations 
lat = c(54.23, 55.14, 55.14, 55.59, 55.49, 55.38, 55.15, 55.14, 55.25, 55.17, 55.16, 54.86,
54.80, 54.95, 54.82, 54.80, 54.80, 54.77, 54.76, 55.48, 55.48, 54.56, 54.55, 54.54, 54.50,
54.63, 54.59, 54.52, 54.40, 54.37, 54.36, 54.16, 55.47, 55.46, 55.12, 55.43, 55.52, 55.62,
55.58, 55.47, 55.35, 55.30, 55.33, 55.32, 55.17, 54.63, 54.95, 54.94, 54.71, 54.36)

long = c(2.730, 1.329, 1.329, 3.225, 1.954, 1.833, 2.090, 2.085, 1.956, 1.643, 1.641, 2.089,
         2.336, 1.489, 1.180, 1.493, 1.493, 1.960, 1.958, 2.559, 2.559, 1.344, 1.343, 1.498,
         1.652, 2.090, 2.331, 2.089, 1.844, 2.335, 2.335, 2.084, 2.903, 2.904, 2.335, 2.335,
         2.338, 2.340, 1.949, 1.469, 1.483, 1.484, 2.901, 2.901, 2.897, 1.040, 1.024, 2.738,
         2.737, 2.551)

nspecies = c(28,16,22,23,17,13,28,18,20,41,21,14,19,41,28,4,32,31,16,9,14,6,35,
          18,9,35,23,5,18,27,27,16,22,16,
          29,11,8,23,28,23,18,16,16,47,31,17,13,23,19,20)

u = c(0,0.1,0.3,0.5,0.7,1,1.5,2.4)

semiv = svariog(long, lat, nspecies, u)

par(mfrow=c(2,2))
plot(semiv$dcloud, semiv$zcloud, xlab='Distance', ylab='Cloud')
plot(semiv$mid, semiv$cla, xlab='Distance', ylab='Classical')
plot(semiv$mid, semiv$med, xlab='Distance', ylab='Median')
plot(semiv$mid, semiv$rob, xlab='Distance', ylab='Robust')
}
\keyword{semi-variogram}
