\name{CorrectSamplingEffort}
\alias{CorrectSamplingEffort}

\title{
Corrects abundance measures for sampling effort
}
\description{
Corrects any measures of abundance that are effort-sensitive according to recorded sampling effort
}
\usage{
CorrectSamplingEffort(diversity)
}

\arguments{
  \item{diversity}{
The PREDICTS species-level data-frame, or a subset thereof
}
}
\details{

}
\value{
Returns the PREDICTS data-frame, but where the column 'Measurement' is corrected for sampling effort.
}
\references{
Newbold et al. (2014). A global model of the response of tropical and sub-tropical forest biodiversity to anthropogenic pressures. <i>Proceedings of the Royal Society, Series B</i> <b>281</b>: 20141371.
}
\author{
Tim Newbold <t.newbold@ucl.ac.uk>
}
\note{

}

\seealso{

}
\examples{
predicts <- CorrectSamplingEffort(predicts)
}
