\name{GVFC}

\alias{GVFC}

\title{
Calculates the raw Visual Fast Count (VFC) estimator of the mean abundance per transect
}
\description{
The function considers the counts per segment and uses them sequentially until \code{d} positive counts
are obtained or until all \code{s} segments have been considered. If we assume that \code{u} counts are used (of
which some may be zero) then
the visual fast count estimator is the mean over the \code{u} counts multiplied by \code{s}. This function is used
by \code{\link{GVFCMOM}} to obtain the method of moments VFC estimator - which has reduced bias compared to the
raw VFC estimator.
}

\usage{
GVFC(counts, s, d) 
}

\arguments{
  \item{counts}{ Vector of length s that contains a count for each segment }
  \item{s}{ Number of segments }
  \item{d}{ Number of positive segment counts needed before counting stops }
}
\value{The raw VFC estimate of the segment mean}

\references{
Barry J and Coggan R (2010) The Visual Fast Count method: a critical examination and development for
underwater video sampling. Aquatic Biology, 11, 101-112.
}
\author{
Jon Barry: jon.barry@cefas.co.uk
}
\note{
The reference given assumes a value of d=1
}

\seealso{\code{\link{GVFCMOM}} }