MergeSites <-
function(diversity, silent=FALSE, match.extra, merge.extra) {
  # Sites that share the same values of match.cols are canididates for merging
  
  match.cols <- c('Source_ID','Reference',
                  'Study_number','Study_name',
                  'Diversity_metric', 
                  'Diversity_metric_unit',
                  'Diversity_metric_type',
                  'Diversity_metric_is_effort_sensitive',
                  'Diversity_metric_is_suitable_for_Chao',
                  'Sampling_method','Sampling_effort_unit',
                  'Study_common_taxon','Rank_of_study_common_taxon',
                  'Longitude','Latitude','Predominant_land_use', 
                  'Use_intensity', 'Sample_start_earliest', 'Sample_end_latest',
                  'Sample_date_resolution', 'Block')
  
  if(!missing(match.extra)) match.cols <- c(match.cols, match.extra)
  
  # Sites to be merged may have different values of site.level.merge.cols - 
  # a message is printed for each column for which there is more than one value 
  # among sites. If 'ignore.duplicates' is TRUE, the first of the non-unique 
  # values is used in the returned data.frame. If 'ignore.duplicates' is FALSE 
  # then the function raises an error once all candidate sites have been 
  # examined.
  
  site.level.merge.cols <- c('Country_distance_metres','Coordinates_method',
                             'Max_linear_extent_metres','Habitat_patch_area_square_metres',
                             'Sampling_effort','Rescaled_sampling_effort',
                             'Habitat_as_described',
                             'Source_for_predominant_land_use',
                             'Km_to_nearest_edge_of_habitat',
                             'Years_since_fragmentation_or_conversion',
                             'Transect_details','Sample_midpoint',
                             'Ecoregion_distance_metres',
                             'Country','UN_subregion','UN_region',
                             'Ecoregion','Wilderness_area',
                             'Biome','Realm','Hotspot')
  
  if(!missing(merge.extra)) site.level.merge.cols <- c(site.level.merge.cols, merge.extra)
  # TODO Take mean?
  
  # Taxa within sites to be merged must have the same values of 
  # taxon.level.unique.cols
  
  taxon.level.unique.cols <- c('Taxon_number','Taxon_name_entered',
                               'Indication','Parsed_name',
                               'Best_guess_binomial','COL_ID','Taxon','Name_status',
                               'Rank','Kingdom','Phylum','Class','Order','Family','Genus','Species',
                               'Higher_taxon')
  
  # Values of synthesize.cols are recomputed after sites have been merged
  
  synthesize.cols <- c('SS','Site_number','Site_name',
                       'SSS','SSB','SSBS',
                       'Measurement','Effort_corrected_measurement')
  
  stopifnot(all.equal(synthesize.cols,setdiff(colnames(diversity), c(match.cols, site.level.merge.cols, taxon.level.unique.cols))))
  
  # Identify studies with duplicated sites
  # Same coordinates, sampling dates, habitat, use intensity
  # Exact match of coordinates for now. In the future, use coordinates precision
  sites <- diversity[!duplicated(diversity$SSS),c(synthesize.cols,match.cols)]
  sites.within.studies <- split(sites[,match.cols], sites$SS)
  candidates <- sapply(sites.within.studies, function(rows) nrow(rows)!=nrow(unique(rows)))
  candidates <- names(which(candidates))
  
  non.unique.site.level.values <- FALSE
  rows.to.drop <- NULL
  
  # Merge
  for(study in candidates) {
    if (!silent) cat('Examining [', study, ']\n', sep='')
    rows <- sites[sites$SS==study,]
    dm.type <- as.character(rows$Diversity_metric_type[1])
    
    # Paste together those columns that should be unique
    ids <- sapply(1:nrow(rows), function(n) paste(do.call('c', rows[n,match.cols]), collapse=' '))
    
    # Split site numbers by id
    consider.sites <- split(sites[sites$SS==study,'Site_number'], ids)
    
    # Consider only sites that appear to be duplicated
    consider.sites <- consider.sites[sapply(consider.sites, length)>1]
    
    # Examine each set of duplicated sites
    for(n in consider.sites) {
      # The rows to be merged
      consider.rows <- which(diversity$SS==study & diversity$Site_number %in% n)
      
      # match.cols must be the same among sites
      stopifnot(1==nrow(unique(diversity[consider.rows,match.cols])))
      
      # Values in site.level.merge.cols are expected to be different among 
      # sites. If there are differences, take values from the first of the 
      # sites to be merged
      x <- unique(diversity[consider.rows,site.level.merge.cols])
      if(nrow(x)>1) {
        contain.duplicates <- sapply(site.level.merge.cols, function(col) 1<length(unique(x[,col])))
        if (!silent){
          cat(paste('Taking the first of non-unique values of [', 
                     names(which(contain.duplicates)), ']', sep=''), sep='\n')}
        non.unique.site.level.values <- TRUE
      }
      
      if (!silent){
        cat('Merging [', study, '] sites [', paste(n, collapse=','), '] with [', 
             length(consider.rows), '] measurements\n', sep='')}
      
      # Values in taxon.level.unique.cols must be unique within taxa within 
      # studies
      stopifnot(all(by(diversity[consider.rows,taxon.level.unique.cols], 
                       diversity[consider.rows,'Taxon_number'], 
                       function(rows) 1==nrow(unique(rows)))))
      
      # Merge measurements
      if(dm.type %in% c('Abundance','Species richness')) {
        # Where Diversity_metric_type is Abundance or Species richness, compute 
        # weighted-mean where weights are sampling efforts
        if (!silent) cat('Computing weighted-mean of', dm.type, '\n')
        measurements <- mapply(FUN=weighted.mean, 
                               split(diversity[consider.rows,'Measurement'],     diversity[consider.rows,'Taxon_number']), 
                               split(diversity[consider.rows,'Sampling_effort'], diversity[consider.rows,'Taxon_number']))
      } else if('Occurrence'==dm.type) {
        # Where Diversity_metric_type is occurrence compute logical 'or'
        if (!silent) cat('Computing logical OR of', dm.type, '\n')
        measurements <- tapply(diversity[consider.rows,'Measurement'], 
                               diversity[consider.rows,'Taxon_number'], 
                               function(m) any(m>0))
        measurements <- ifelse(measurements, 1, 0)
      } else {
        stop(dm.type)
      }
      
      n.taxa <- length(unique(diversity[consider.rows,'Taxon_number']))
      if(sum(diversity[consider.rows,'Site_number']==n[1]) != n.taxa) {
        # Each taxon is not represented at the first site in n - synthesize new 
        # rows.
        if (!silent){
          cat('Not all taxa are represented at the first site. Synthesizing', 
               length(measurements), 'rows\n')}
        
        # How could this ever happen?
        stopifnot(n.taxa<length(consider.rows))
        
        # The taxa that where measured by this study
        taxa <- unique(diversity[consider.rows,taxon.level.unique.cols])
        taxa <- taxa[order(taxa$Taxon_number),]
        
        # Set the first n.taxa rows to the contents of the first site
        diversity[consider.rows[1:n.taxa],] <- diversity[consider.rows[1],]
        
        # Set the taxon columns for the first n.taxa rows
        diversity[consider.rows[1:n.taxa],taxon.level.unique.cols] <- taxa
      }
      
      if (!silent) cat('Replacing', length(measurements), 'merged measurements\n')
      overwrite <- consider.rows[diversity[consider.rows,'Site_number']==n[1]]
      stopifnot(length(overwrite)==length(measurements))
      diversity[overwrite,'Measurement'] <- unname(measurements)
      
      # Values have been merged to the first of the site numbers given in n. 
      # Drop the other sites. 
      drop <- consider.rows[diversity[consider.rows,'Site_number']!=n[1]]
      stopifnot(length(drop)+length(overwrite)==length(consider.rows))
      rows.to.drop <- c(rows.to.drop, drop)
    }
  }
  
  if(!is.null(rows.to.drop)) {
    if (!silent) cat('Dropping', length(rows.to.drop), 'measurements that have been merged\n')
    diversity <- diversity[-rows.to.drop,]
    
    if (!silent) cat('Dropping unused factor levels\n')
    diversity <- droplevels(diversity)
  }
  
  return (diversity)
  
}
