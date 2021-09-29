\name{SiteMetrics}
\alias{SiteMetrics}

\title{
Calculates site-level biodiversity metrics
}
\description{
Calculates species richness and total abundance at each PREDICTS sites
}
\usage{
SiteMetrics(diversity, extra.cols = NULL,srEstimators = c("Chao","Rare"))
}

\arguments{
  \item{diversity}{
The raw species-level PREDICTS diversity data
}
  \item{extra.cols}{
Any non-standard columns to be transferred to the final site-level data frame
}
  \item{srEstimators}{
Vector of species richness estimators to calculate. Recognizes "Chao" for the Chao species richness estimator (Chao et al., 2005) or "Rare" for rarefied richness.  
}
}
\details{

}
\value{
A data frame where each row corresponds to a PREDICTS site
}
\references{
Chao A., Chazdon R.L., Colwell R.K. & Shen T.-J. (2005). A new statistical approach for assessing similarity of species composition with incidence and abundance data. <i>Ecology Letters</i> <b>8</b>: 148-159.
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

