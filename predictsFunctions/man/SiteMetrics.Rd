\name{SiteMetrics}
\alias{SiteMetrics}

\title{
Calculates site-level biodiversity metrics
}
\description{
Calculates species richness and total abundance at each PREDICTS sites
}
\usage{
SiteMetrics(diversity, extra.cols = NULL)
}

\arguments{
  \item{diversity}{
The raw species-level PREDICTS diversity data
}
  \item{extra.cols}{
Any non-standard columns to be transferred to the final site-level data frame
}
}
\details{

}
\value{
A data frame where each row corresponds to a PREDICTS site
}
\references{
Hudson L.N., Newbold, T. et al. (2014). The PREDICTS database: a global database of how local terrestrial biodiversity responds to human impacts. <i>Ecology & Evolution</i> <b>4</b>: 4701-4735.
}
\author{
Tim Newbold <t.newbold@ucl.ac.uk>
}
\note{

}

\seealso{

}
\examples{
sites <- SiteMetrics(diversity = predicts)
}

