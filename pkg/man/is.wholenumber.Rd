\name{is.wholenumber}
\alias{is.wholenumber}

\title{
To check whether an argument is an integer
}

\description{
Used in error checking to ascertain whether a function argument is an integer.
}

\usage{
is.wholenumber(x, tol = .Machine$double.eps^0.5)
}

\arguments{
\item{x}{Number to be checked.}
\item{tol}{If \code{x} is closer to an integer than this, then it passes.}
}

\value{
Vector of logical values if the corresponding input values is an integer or not.
}

\references{
\code{is.wholenumber} is taken from the \code{\link[base]{is.integer}} help file.
}

\examples{
is.wholenumber( seq(1, 5, by = 0.5) ) #-->  TRUE FALSE TRUE ...
}
