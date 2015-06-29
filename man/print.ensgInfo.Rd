\name{print.ensgInfo}
\docType{methods}
\alias{print,ensgInfo-method}
\alias{print.ensgInfo}

\title{Print an ensgInfo Object}
\description{
  Prints an \code{ensgInfo} object.
}
\usage{
 \method{print}{ensgInfo}(x, full=FALSE, ...)
}

\arguments{
  \item{x}{Object of class \code{ensgInfo}.}
  \item{full}{Logical, shall the full information be plotted.}
  \item{...}{Additional arguments.}
}

\details{
The print function displays an \code{ensgInfo} object. By default just the Ensembl ID and the corresponding gene name is
plotted. Setting the option \code{full=TRUE} provides further information.
}

\author{Daniel Fischer}

\keyword{methods}
\keyword{print}
