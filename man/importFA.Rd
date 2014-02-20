\name{importFA}
\alias{importFA}

\title{
 Importing a Fasta File.
}

\description{
  This function imports a standard fasta file.
}
\usage{
  importFA(file)
}

\arguments{
  \item{file}{Specifies the filename/path.}
}

\details{
 This function imports a standard fasta file. It assumes that label and sequence lines are alternating,
 meaning in the odd lines is the sequence name given, starting with > and in the even rows are the 
 corresponding sequences.
}
\value{
  A vector containg the sequences. The vector names coorespond to the sequence names given in the fasta file.
}

\author{
  Daniel Fischer
}

\examples{
\dontrun{
importFA("/home/data/myFasta.fa")
}
}

\keyword{methods}
 
