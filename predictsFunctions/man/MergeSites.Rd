\name{MergeSites}
\alias{MergeSites}

\title{
Merge PREDICTS sites having the same coordinates
}
\description{
Merges PREDICTS sites that have the same coordinates, given that these generally reflect multiple sub-samples within one site (e.g. many pit-fall traps along a single transect)
}
\usage{
MergeSites(diversity, silent = FALSE, match.extra, merge.extra)
}

\arguments{
  \item{diversity}{
The PREDICTS species-level data-frame, or a subset thereof
}
  \item{silent}{
Suppress progress reporting
}
  \item{match.extra}{
Extra columns that are required to match for merging (see below for details)
}
  \item{merge.extra}{
Extra columns to be merged (see below for details)
}
}
\details{
If values in some columns match ('match' columns), sites are candidates for merging. These columns are: 'Source_ID', 'Reference', 'Study_number', 'Study_name', 'Diversity_metric, 'Diversity_metric_unit', 'Diversity_metric_type', 'Diversity_metric_is_effort_sensitive, 'Diversity_metric_is_suitable_for_Chao', 'Sampling_method', 'Sampling_effort_unit', 'Study_common_taxon', 'Rank_of_study_common_taxon', 'Longitude', 'Latitude', 'Predominant_land_use', 'Use_intensity', 'Sample_start_earliest', 'Sample_end_latest','Sample_date_resolution', and 'Block'.

Other 'merge' columns do not have to match for merging to occur, and the value from the first site is taken upon merging of two sites. Default merge columns are 'Country_distance_metres', 'Coordinates_method', 'Max_linear_extent_metres', 'Habitat_patch_area_square_metres', 'Sampling_effort', 'Rescaled_sampling_effort', 'Habitat_as_described', 'Source_for_predominant_land_use', 'Km_to_nearest_edge_of_habitat', 'Years_since_fragmentation_or_conversion', 'Transect_details', 'Sample_midpoint', 'Ecoregion_distance_metres', 'Country', 'UN_subregion', 'UN_region', 'Ecoregion', 'Wilderness_area', 'Biome', 'Realm', and 'Hotspot'

Additional 'match' and 'merge' columns can be specified using the 'match.extra' and 'merge.extra' parameters.
}
\value{
Returns a copy of the PREDICTS data-frame where qualifying sites have been merged together
}
\references{
Newbold et al. (2015). Global effects of land use on local terrestrial biodiversity. <i>Nature</i> <b>520</b>: 45-50.
}
\author{
Tim Newbold <t.newbold@ucl.ac.uk>
}
\note{

}

\seealso{

}
\examples{
predicts <- MergeSite(predicts)
}
