\name{importGTF}
\alias{importGTF}

\title{
  Import a GTF File
}

\description{
  This function imports a gtf file.
}
\usage{
  importGTF(gtf)
}

\arguments{
  \item{gtf}{File name of the gtf file}
}

\details{
  This function imports a gtf file and splits the column 9 which is usually tricky to handle as the order of the variables is not always
  the same.
}
\value{
  A data frame containing the columns of the gtf file, including the splitted last column.
}

\author{
  Daniel Fischer
}

 
