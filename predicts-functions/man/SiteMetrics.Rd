\name{SiteMetrics}
\alias{SiteMetrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SiteMetrics(diversity, extra.cols = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{diversity}{
%%     ~~Describe \code{diversity} here~~
}
  \item{extra.cols}{
%%     ~~Describe \code{extra.cols} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (diversity, extra.cols = NULL) 
{
    cat(paste("Computing site metrics for", nrow(diversity), 
        "measurements\n"))
    stopifnot(all(diversity$Diversity_metric_is_valid))
    stopifnot(all(!is.na(diversity$Sampling_effort)))
    site.cols <- c("Source_ID", "Study_number", "Study_name", 
        "Site_number", "Site_name", "Block", "Predominant_habitat", 
        "Use_intensity", "Longitude", "Latitude", "Sample_start_earliest", 
        "Sample_end_latest", "Diversity_metric_type")
    if (!is.null(extra.cols)) 
        site.cols <- union(site.cols, extra.cols)
    site.cols <- intersect(site.cols, colnames(diversity))
    cat("The data contain", nlevels(diversity$Source_ID), "sources,", 
        nlevels(diversity$SS), "studies and", nlevels(diversity$SSS), 
        "sites\n")
    if (any("Diversity index" == diversity$Diversity_metric_type)) {
        stop("Diversity indices are not yet supported")
    }
    bad <- setdiff(diversity$Diversity_metric_type, c("Abundance", 
        "Occurrence", "Species richness"))
    if (length(bad) > 0) {
        stop("Unrecognied diversity metrics ", paste(bad, collapse = ","))
    }
    cat("Computing site-level values\n")
    diversity$Is_abundance <- "Abundance" == diversity$Diversity_metric_type
    diversity$Is_occurrence <- "Occurrence" == diversity$Diversity_metric_type
    diversity$Is_species_richness <- "Species richness" == diversity$Diversity_metric_type
    site.abundance <- tapply(diversity$Is_abundance, diversity$SSS, 
        unique)
    site.occurrence <- tapply(diversity$Is_occurrence, diversity$SSS, 
        unique)
    site.species.richness <- tapply(diversity$Is_species_richness, 
        diversity$SSS, unique)
    total.abundance <- SiteTotalAbundance(diversity, site.abundance)
    species.richness <- SiteSpeciesRichness(diversity, site.abundance, 
        site.occurrence, site.species.richness)
    cat("Assembling site-level values\n")
    res <- cbind(diversity[!duplicated(diversity$SSS), c(site.cols, 
        "SS", "SSS")], Total_abundance = total.abundance, Species_richness = species.richness)
    rownames(res) <- NULL
    res <- res[order(res$Source_ID, res$Study_number, res$Site_number), 
        ]
    stopifnot(all.equal(as.character(res$SS), paste(res$Source_ID, 
        res$Study_number)))
    stopifnot(all.equal(as.character(res$SSS), paste(res$Source_ID, 
        res$Study_number, res$Site_number)))
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
