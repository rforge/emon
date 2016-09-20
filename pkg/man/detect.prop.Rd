\name{detect.prop}
\alias{detect.prop}

\title{
Probability of detecting a feature that covers a proportion \code{theta} of the survey area.
}
\description{
The function can calculate the probability of a feature that occupies
a proportion \code{theta} of the sampling area and where the sampling point density of the survey is
specified; the sampling point density needed to achieve a specified probability of detection, where
\code{theta} is also specified ; or the value of \code{theta} that will be detected with specified
probability and sampling density.Unless the feature is made of a large number of random segments
(see below for how to deal with this situation), these methods apply only when the pattern of points
in the sampling deisgn is random.
}
\usage{
detect.prop(statistic, theta=NA, pdetect=NA, ssize=NA)
}
\arguments{
  \item{statistic}{
Describes what aspect of design you want calculated. The choices are \code{"P"} (probability detection); 
\code{"N"} (sample size) or \code{"F"} (feature proportion).
}
  \item{theta}{
Feature proportion. Not needed if \code{statistic="F"}.
}
  \item{pdetect}{
Probability detection. Not needed if \code{statistic="P"}.
}
  \item{ssize}{
Sample size. Not needed if \code{statistic="N"}.
}
}

\details{ The probability of detection is \code{p = 1 - (1 - theta)^{N}}. Formulae for \code{theta} and
\code{N} are readily obtained from this formula. If the spatial pattern of the feature consists of lots of small,
random uniformly distributed fragments, then we can redefine \code{theta = Na/A} where \code{a} is the
area of the sampling unit and \code{A} is the sampling area.In this situation, the probability of patch detection
applies no matter what the spatial pattern of points in the sampling design. Unlike \code{detect}, \code{detect.prop}
works for vectors - so long as the input vectors are of the same length.
}

\value{
\item{prob}{Probability of detection}
\item{ssize}{Sample size}
\item{prop}{Feature proportion}
}

\author{
Jon Barry: Jon.Barry@cefas.co.uk
}

\examples{
detect.prop(statistic='P', theta=0.02, ssize=80)
detect.prop(statistic='N', theta=0.02, pdetect=0.9)
detect.prop(statistic='F', pdetect=0.9, ssize=80)
}