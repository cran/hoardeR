\name{speciesFigure}
\alias{speciesFigure}

\title{
 Showing Quantities of Different Species.
}

\description{
  This function visualizes the different quantities of blast machtes.
}
\usage{
  speciesFigure(xml, species=NULL, type="chr", n=2:11, plot=TRUE)
}

\arguments{
  \item{xml}{An xml object.}
  \item{species}{A vector with species names.}
  \item{type}{The type of plot.}
  \item{n}{Ranks to be plotted}
  \item{plot}{Logical, shall the plot be plotted}
}

\details{
 This function plots the frequency barplot of the blast results, divided for each species. The species of interest can be
 provided with the \code{species} object. 
}
\value{
  A figure.
}

\author{
  Daniel Fischer
}

\examples{
\dontrun{
  speciesFigure(xml, species=NULL, type="chr", n=2:11, plot=TRUE)
}
}

\keyword{methods}
 
