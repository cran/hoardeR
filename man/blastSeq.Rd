\name{blastSeq}
\alias{blastSeq}

\title{
 Sending a Fasta Sequence to NCBI Blast service
}

\description{
  This function sends fasta sequences to the NCBI Blast service.
}
\usage{
  blastSeq(seq, n_blast=20, delay_req=3, delay_rid=60, email=NULL, 
           xmlFolder=NULL, keepInMemory=TRUE)
}

\arguments{
  \item{seq}{The fasta sequence that should be blasted (\code{String}).}
  \item{n_blast}{Amount of parallel blast requests, in case \code{seq} is a vector.}
  \item{delay_req}{Seconds between the single Blast requests.}
  \item{delay_rid}{Seconds between the single result requests.}
  \item{email}{User email, required information from NCBI (\code{String}).}
  \item{xmlFolder}{Path to the result folder.}
  \item{keepInMemory}{Logical, shall the results be kept in the memory.}
}

\details{
 This function sends fasta sequence to the NCBI blast service. The defaults for the delays are required by NCBI. Also, 
 NCBI asks to provide an email address.
 
 In case the option \code{keepInMemory} is set to TRUE the Blast results will be kept in memory, otherwise they will
 be just written to the HDD. 
 
}
\value{
  An xml file that contains the the NCBI result.
}

\author{
  Daniel Fischer
}

\examples{
\dontrun{
blastSeq("ACGTGCATCGACTAGCTACGACTACGACTATC")
}
}

\keyword{methods}
 
