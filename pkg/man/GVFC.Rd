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
Barry J, Eggleton J, Ware S and Curtis M (2015) Generalizing Visual Fast Count Estimators for Underwater
Video Surveys. Ecosphere. 
http://www.esajournals.org/doi/full/10.1890/ES15-00093.1
}
\author{
Jon Barry: jon.barry@cefas.co.uk
}

\seealso{\code{\link{GVFCMOM}} }