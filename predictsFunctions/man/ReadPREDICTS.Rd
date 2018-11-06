\name{ReadPREDICTS}
\alias{ReadPREDICTS}

\title{
Read the PREDICTS database
}
\description{
Reads the PREDICTS database from an .rds file
\usage{
ReadPREDICTS(path)
}

\arguments{
  \item{path}{
The full path to a local copy of the PREDICTS database, which can be downloaded from: http://dx.doi.org/10.5519/0066354.
}
}
\details{

}
\value{
Returns the PREDICTS database as an R data-frame
}
\references{
Hudson L.N., Newbold, T. et al. (2017). The database of the PREDICTS (Projecting Responses of Ecological Diversity In Changing Terrestrial Systems) project. <i>Ecology & Evolution</i> <b>7</b>: 145-188.
}
\author{
Tim Newbold <t.newbold@ucl.ac.uk>
}
\note{

}

\seealso{

}
\examples{
predicts.path <- "C:/Users/tim_n/Documents/database.rds"

predicts <- ReadPREDICTS(predicts.path)
}
