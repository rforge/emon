\name{detect}
\alias{detect}

\title{
Probability of circular patch detection
}
\description{
The function can calculate the probability of detection of a circular patch of specified radius for a 
specified density of points; the density needed to achieve a specified probability of detection; or the
radius of the patch that will be detected with specified probability and sampling density.This is done for
random, square lattice, and triangular lattice spatial sampling designs.
}
\usage{
detect(method, statistic, area=NA, radius=NA, pdetect=NA, ssize=NA)
}
\arguments{
  \item{method}{
Defines the spatial sampling design to be used. The values can be \code{"R"} (random), \code{"S"} (square lattice)
or \code{"T"} (triangular lattice). See Barry and Nicholson (1993) for details and formulae for the 
probabilities of detection for the square
lattice and triangular lattice designs. For the random design, \code{prob(detect)=1 - (1 - a/A)^N}, where
\code{a} is the patch area and \code{A} is the survey area. This gives similar answers to the formula in
Barry and Nicholson, but is exact for fixed sample size.
}
  \item{statistic}{
Describes what aspect of design you want calculated. The choices are \code{"P"} (probability detection); 
\code{"N"} (sample size) or \code{"R"} (patch radius).
}
  \item{area}{
The survey area (same units as distance and radius).
}
  \item{radius}{
Patch radius. Not needed if \code{statistic="R"}.
}
  \item{pdetect}{
Probability detection. Not needed if \code{statistic="P"}.
}
  \item{ssize}{
Sample size. Not needed if \code{statistic="N"}.
}
}

\details{ The basic idea is that you wish to conduct a survey in an area \code{area} to detect some object (patch) of
interest. This could be a cockle patch, an area of reef or an archaeological deposit. This function asssumes that
the object is circular with radius \code{radius}. You have three choices of sampling deign to use: spatial, square
lattice and triangular lattice. In terms of patch detection, for a given sample size, the triangular design gives
the highest probability - because its points are equi-distant apart.

The simplest application of this function is to assess the patch detection probability for a particular design. This
is obtained using the \code{statistic="P"} option. However, the problem can be turned around and this function used to
calculate the sample size needed to obatain a specific patch detection probability (\code{statistic="N"}) or the radius
of the patch that would be detected with some desired probability (\code{statistic="R"). Th}is last scenario might be
useful if there was some particular size of patch that you wanted to be sure (say, 90 percent) of detecting.
}

\value{
\item{prob}{Probability of patch detection}
\item{ssize}{Sample size}
\item{rad}{Patch radius}
\item{sep}{Separation distance (for square and triangular lattice designs)}
}

\references{
Barry J and Nicholson M D (1993) Measuring the probabilities of patch detection
for four spatial sampling schemes. Journal of Applied Statistics, 20, 353-362.                    
}
\author{
Jon Barry: Jon.Barry@cefas.co.uk
}

\examples{
detect(method='R', statistic='P', area=100, radius=2, ssize=15)$prob
detect(method='R', statistic='N', area=100, radius=2, pdetect=0.95)$ssize
detect(method='R', statistic='R', area=100, pdetect=0.95, ssize=15)$rad

detect(method='S', statistic='P', area=100, radius=1.4, ssize=15)
detect(method='S', statistic='N', area=100, radius=1.4, pdetect=0.6)

# Plot patch detection as a function of radius
square = rep(0,200); rand = square; triang = rand
radius = seq(0.01, 2, 0.01)

for (j in 1:200) {
rand[j] = detect(method='R', statistic='P', area=100, radius=radius[j], ssize=15)$p
square[j] = detect(method='S', statistic='P', area=100, radius=radius[j], ssize=15)$p
triang[j] = detect(method='T', statistic='P', area=100, radius=radius[j], ssize=15)$p
}

plot(radius, rand, ylim=c(0,1), xlab='Patch radius', ylab='Probability detection', type='l')
lines(radius, square, col=2, lty=2)
lines(radius, triang, col=3, lty=3)
legend('topleft', legend=c('Random', 'Square', 'Triangular'), col=c(1,2,3), lty=c(1,2,3))
}