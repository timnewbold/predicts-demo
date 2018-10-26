\name{MergeSites}
\alias{MergeSites}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MergeSites(diversity, silent = FALSE, match.extra, merge.extra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{diversity}{
%%     ~~Describe \code{diversity} here~~
}
  \item{silent}{
%%     ~~Describe \code{silent} here~~
}
  \item{match.extra}{
%%     ~~Describe \code{match.extra} here~~
}
  \item{merge.extra}{
%%     ~~Describe \code{merge.extra} here~~
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
function (diversity, silent = FALSE, match.extra, merge.extra) 
{
    match.cols <- c("Source_ID", "Reference", "Study_number", 
        "Study_name", "Diversity_metric", "Diversity_metric_unit", 
        "Diversity_metric_type", "Diversity_metric_is_effort_sensitive", 
        "Diversity_metric_is_suitable_for_Chao", "Sampling_method", 
        "Sampling_effort_unit", "Study_common_taxon", "Rank_of_study_common_taxon", 
        "Longitude", "Latitude", "Predominant_land_use", "Use_intensity", 
        "Sample_start_earliest", "Sample_end_latest", "Sample_date_resolution", 
        "Block")
    if (!missing(match.extra)) 
        match.cols <- c(match.cols, match.extra)
    site.level.merge.cols <- c("Country_distance_metres", "Coordinates_method", 
        "Max_linear_extent_metres", "Habitat_patch_area_square_metres", 
        "Sampling_effort", "Rescaled_sampling_effort", "Habitat_as_described", 
        "Source_for_predominant_land_use", "Km_to_nearest_edge_of_habitat", 
        "Years_since_fragmentation_or_conversion", "Transect_details", 
        "Sample_midpoint", "Ecoregion_distance_metres", "Country", 
        "UN_subregion", "UN_region", "Ecoregion", "Wilderness_area", 
        "Biome", "Realm", "Hotspot")
    if (!missing(merge.extra)) 
        site.level.merge.cols <- c(site.level.merge.cols, merge.extra)
    taxon.level.unique.cols <- c("Taxon_number", "Taxon_name_entered", 
        "Indication", "Parsed_name", "Best_guess_binomial", "COL_ID", 
        "Taxon", "Name_status", "Rank", "Kingdom", "Phylum", 
        "Class", "Order", "Family", "Genus", "Species", "Higher_taxon")
    synthesize.cols <- c("SS", "Site_number", "Site_name", "SSS", 
        "SSB", "SSBS", "Measurement", "Effort_corrected_measurement")
    stopifnot(all.equal(synthesize.cols, setdiff(colnames(diversity), 
        c(match.cols, site.level.merge.cols, taxon.level.unique.cols))))
    sites <- diversity[!duplicated(diversity$SSS), c(synthesize.cols, 
        match.cols)]
    sites.within.studies <- split(sites[, match.cols], sites$SS)
    candidates <- sapply(sites.within.studies, function(rows) nrow(rows) != 
        nrow(unique(rows)))
    candidates <- names(which(candidates))
    non.unique.site.level.values <- FALSE
    rows.to.drop <- NULL
    for (study in candidates) {
        if (!silent) 
            cat("Examining [", study, "]\n", sep = "")
        rows <- sites[sites$SS == study, ]
        dm.type <- as.character(rows$Diversity_metric_type[1])
        ids <- sapply(1:nrow(rows), function(n) paste(do.call("c", 
            rows[n, match.cols]), collapse = " "))
        consider.sites <- split(sites[sites$SS == study, "Site_number"], 
            ids)
        consider.sites <- consider.sites[sapply(consider.sites, 
            length) > 1]
        for (n in consider.sites) {
            consider.rows <- which(diversity$SS == study & diversity$Site_number \%in\% 
                n)
            stopifnot(1 == nrow(unique(diversity[consider.rows, 
                match.cols])))
            x <- unique(diversity[consider.rows, site.level.merge.cols])
            if (nrow(x) > 1) {
                contain.duplicates <- sapply(site.level.merge.cols, 
                  function(col) 1 < length(unique(x[, col])))
                if (!silent) {
                  cat(paste("Taking the first of non-unique values of [", 
                    names(which(contain.duplicates)), "]", sep = ""), 
                    sep = "\n")
                }
                non.unique.site.level.values <- TRUE
            }
            if (!silent) {
                cat("Merging [", study, "] sites [", paste(n, 
                  collapse = ","), "] with [", length(consider.rows), 
                  "] measurements\n", sep = "")
            }
            stopifnot(all(by(diversity[consider.rows, taxon.level.unique.cols], 
                diversity[consider.rows, "Taxon_number"], function(rows) 1 == 
                  nrow(unique(rows)))))
            if (dm.type \%in\% c("Abundance", "Species richness")) {
                if (!silent) 
                  cat("Computing weighted-mean of", dm.type, 
                    "\n")
                measurements <- mapply(FUN = weighted.mean, split(diversity[consider.rows, 
                  "Measurement"], diversity[consider.rows, "Taxon_number"]), 
                  split(diversity[consider.rows, "Sampling_effort"], 
                    diversity[consider.rows, "Taxon_number"]))
            }
            else if ("Occurrence" == dm.type) {
                if (!silent) 
                  cat("Computing logical OR of", dm.type, "\n")
                measurements <- tapply(diversity[consider.rows, 
                  "Measurement"], diversity[consider.rows, "Taxon_number"], 
                  function(m) any(m > 0))
                measurements <- ifelse(measurements, 1, 0)
            }
            else {
                stop(dm.type)
            }
            n.taxa <- length(unique(diversity[consider.rows, 
                "Taxon_number"]))
            if (sum(diversity[consider.rows, "Site_number"] == 
                n[1]) != n.taxa) {
                if (!silent) {
                  cat("Not all taxa are represented at the first site. Synthesizing", 
                    length(measurements), "rows\n")
                }
                stopifnot(n.taxa < length(consider.rows))
                taxa <- unique(diversity[consider.rows, taxon.level.unique.cols])
                taxa <- taxa[order(taxa$Taxon_number), ]
                diversity[consider.rows[1:n.taxa], ] <- diversity[consider.rows[1], 
                  ]
                diversity[consider.rows[1:n.taxa], taxon.level.unique.cols] <- taxa
            }
            if (!silent) 
                cat("Replacing", length(measurements), "merged measurements\n")
            overwrite <- consider.rows[diversity[consider.rows, 
                "Site_number"] == n[1]]
            stopifnot(length(overwrite) == length(measurements))
            diversity[overwrite, "Measurement"] <- unname(measurements)
            drop <- consider.rows[diversity[consider.rows, "Site_number"] != 
                n[1]]
            stopifnot(length(drop) + length(overwrite) == length(consider.rows))
            rows.to.drop <- c(rows.to.drop, drop)
        }
    }
    if (!is.null(rows.to.drop)) {
        if (!silent) 
            cat("Dropping", length(rows.to.drop), "measurements that have been merged\n")
        diversity <- diversity[-rows.to.drop, ]
        if (!silent) 
            cat("Dropping unused factor levels\n")
        diversity <- droplevels(diversity)
    }
    return(diversity)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
