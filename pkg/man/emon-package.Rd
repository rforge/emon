\name{emon-package}

\alias{emon-package}

\alias{emon}

\docType{package}

\title{Tools for environmental and ecological survey design and analysis}

\description{
This package gives six tools for designing and analysing ecological and environmental surveys. The tools
are mainly
designed for marine  and benthic ecology applications, but they could easily be adopted for terrestrial
ecology. Three
of the tools give statistical power for specific survey designs (\code{\link{power.BACI}},
\code{\link{power.groups}} and
\code{\link{power.trend}}). The other three are for more specialised applications. These are: the
generalised visual fast
count estimator for underwater video surveys (\code{\link{GVFCMOM}}); an estimate of the empirical
semi-variogram for
examining spatial correlation between stations (\code{\link{svariog}}); and the patch detection
probability for three
spatial sampling designs (\code{\link{detect}}).
}
\details{
  \tabular{ll}{
    Package: \tab emon\cr
    Type: \tab Package\cr
    Version: \tab 0.6\cr
    Date: \tab 2015-03-31\cr
    License: \tab GPL-3\cr
  }
The six tools in this package are as follows:

Power for BACI designs (\code{\link{power.BACI}}, \code{\link{generate.trend}}, \code{\link{addnoise}},
          \code{\link{mannkendall}}, \code{\link{mannkendall.stat}}, \code{\link{permute.BACI}}).

Power for comparing two groups (\code{\link{power.groups}}, \code{\link{permute.groups}},
          \code{\link{size2.samevar}}).

Power for detecting trends (\code{\link{power.trend}}, \code{\link{generate.trend}}, \code{\link{addnoise}}).

Sample sizes and probabilities for patch detection with different spatial sampling patterns
    (\code{\link{detect}}, \code{\link{fS.detect}}, \code{\link{fT.detect}}).

Semi-variogram function for investigating spatial dependency (\code{\link{svariog}}).

Method of moments estimator for Generalised Visual Fast Count estimation for video surveys
(\code{\link{GVFCMOM}}, \code{\link{GVFC}}, \code{\link{expected.pois}},
 \code{\link{expected.nb}}, \code{\link{mom.min.pois}}, \code{\link{mom.min.nb}}).

The help functions for the individual functions describe the methods used. However, perhaps
the unique feature of
the power functions in \code{emon} is that the statistical power is calculated by simulation.
This has the disadvantage
of increased computing time; however, the advantage is that the power calculations does not rely
on the assumptions
behind many of the theoretical results. The simulation method also means that power can be
calculated for a range of
data distributions and for a variety of statistical tests that might be used to evaluate p-values.
}

\author{
  Jon Barry and David Maxwell
  
  Maintainer: Jon Barry: Jon.Barry@cefas.co.uk
}

\keyword{power, ecological surveys, video surveys, spatial correlation, patch detection}
  
\seealso{ \code{\link{power.BACI}}, \code{\link{power.groups}}, \code{\link{power.trend}},
          \code{\link{detect}}, \code{\link{svariog}}, \code{\link{GVFCMOM}}}

