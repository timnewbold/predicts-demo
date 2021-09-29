SiteMetrics <-
function(diversity, extra.cols=NULL) {
  # Takes a data.frame of diversity measurements and returns a data.frame 
  # of site-level matrics.
  # diversity should have columns
  #   Source_ID
  #   Study_number
  #   Site_number
  #   SS (Source_ID Study_number)
  #   SSS (Source_ID Study_number Site_number)
  #   Diversity_metric_is_valid
  #   Diversity_metric_type
  #   Diversity_metric_is_effort_sensitive
  #   Diversity_metric_is_suitable_for_Chao
  #   Sampling_effort
  #   Sampling_method_is_valid
  #   Measurement
  #   All columns in 'traits'
  
  # The Diverity_* and Sampling_* columns are each guaranteed to be unique 
  # within studies. Because sites are nested within studies, these columns are 
  # also guaranteed to be unique within sites.
  cat(paste('Computing site metrics for', nrow(diversity), 'measurements\n'))
  
  stopifnot(all(diversity$Diversity_metric_is_valid))
  stopifnot(all(!is.na(diversity$Sampling_effort)))
  
  site.cols <- c('Source_ID','Study_number','Study_name','Site_number',
                 'Site_name','Block','Predominant_habitat','Use_intensity',
                 'Longitude','Latitude','Sample_start_earliest',
                 'Sample_end_latest','Diversity_metric_type')
  
  if(!is.null(extra.cols)) site.cols <- union(site.cols, extra.cols)
  
  # Drop names of columns that are not present
  site.cols <- intersect(site.cols, colnames(diversity))
  
  cat('The data contain', nlevels(diversity$Source_ID), 'sources,', 
       nlevels(diversity$SS), 'studies and',
       nlevels(diversity$SSS), 'sites\n')
  
  
  if(any('Diversity index'==diversity$Diversity_metric_type)) {
    stop('Diversity indices are not yet supported')
  }
  
  bad <- setdiff(diversity$Diversity_metric_type, c('Abundance','Occurrence','Species richness'))
  if(length(bad)>0) {
    stop('Unrecognied diversity metrics ', paste(bad, collapse=','))
  }
  
  cat('Computing site-level values\n')
  
  # Measurement-level values
  diversity$Is_abundance <- 'Abundance'==diversity$Diversity_metric_type
  diversity$Is_occurrence <- 'Occurrence'==diversity$Diversity_metric_type
  diversity$Is_species_richness <- 'Species richness'==diversity$Diversity_metric_type
  
  # Site-level values
  site.abundance <- tapply(diversity$Is_abundance, diversity$SSS, unique)
  site.occurrence <- tapply(diversity$Is_occurrence, diversity$SSS, unique)
  site.species.richness <- tapply(diversity$Is_species_richness, diversity$SSS, unique)
  
  total.abundance <- .SiteTotalAbundance(diversity, site.abundance)
  species.richness <- .SiteSpeciesRichness(diversity, site.abundance, 
                                           site.occurrence, 
                                           site.species.richness)
  
  # Calculate specified species richness estimators
  if ("Chao" %in% srEstimators){
    chao <- .SiteChao(diversity, site.abundance)
  }
  
  if ("Rare" %in% srEstimators){
    rsrich<-.SiteRarefiedRichness(diversity,site.abundance)
  }
  
  cat('Assembling site-level values\n')
  res <- cbind(diversity[!duplicated(diversity$SSS),c(site.cols,'SS','SSS')], 
               Total_abundance=total.abundance, 
               Species_richness=species.richness,
               ChaoR=if ("Chao" %in% srEstimators) chao else NA,
               Richness_rarefied=if ("Rare" %in% srEstimators) rsrich else NA)
  rownames(res) <- NULL
  
  res <- res[order(res$Source_ID, res$Study_number, res$Site_number),]
  
  # Sanity checks
  stopifnot(all.equal(as.character(res$SS), 
                      paste(res$Source_ID, res$Study_number)))
  stopifnot(all.equal(as.character(res$SSS), 
                      paste(res$Source_ID, res$Study_number, res$Site_number)))
  
  return (res)
  
}

.SiteTotalAbundance <- function(diversity, site.abundance) {
  cat("Computing total abundance\n")
  ta <- rep(NA, length(site.abundance))
  ta[site.abundance] <- tapply(diversity$Measurement[diversity$Is_abundance], 
                               droplevels(diversity$SSS[diversity$Is_abundance]), 
                               sum)
  return (ta)
}

.SiteSpeciesRichness <- function(diversity, site.abundance, site.occurrence, 
                                site.species.richness) {
  
  interesting <- diversity$Is_abundance | diversity$Is_occurrence
  
  cat('Computing species richness\n')
  # The correct way - taxa are unique within sites
  # The number of non-zero abundances and/or counts
  sr <- rep(NA, length(site.abundance))
  sr[site.abundance | site.occurrence] <- 
    tapply(diversity$Measurement[interesting], 
           droplevels(diversity$SSS[interesting]), 
           function(m) sum(m>0))
  
  sr[site.species.richness] <- 
    tapply(diversity$Measurement[diversity$Is_species_richness], 
           droplevels(diversity$SSS[diversity$Is_species_richness]), 
           sum)
  
  return (sr)
}

.SiteChao <- function(diversity, site.abundance) {
  .Log("Computing Chao\n")
  # Compute Chao only for studies where Diversity_metric_is_suitable_for_Chao 
  # is TRUE and all measurements within the study are integers.
  
  # http://stackoverflow.com/questions/3476782/how-to-check-if-the-number-is-integer
  # TODO all.equal
  study.chao.suitable <- 
    tapply(diversity$Measurement[diversity$Diversity_metric_is_suitable_for_Chao],
           diversity$SS[diversity$Diversity_metric_is_suitable_for_Chao], 
           function(m) all(floor(m)==m))
  
  # Studies for which Diversity_metric_is_suitable_for_Chao is FALSE will have 
  # a value of NA in study.chao.suitable
  study.chao.suitable[is.na(study.chao.suitable)] <- FALSE
  site.chao.suitable <- study.chao.suitable[tapply(diversity$SS, diversity$SSS, unique)]
  
  diversity$Is_Chao_suitable <- site.chao.suitable[diversity$SSS]
  
  # Only abundance metrics should be marked as suitable for Chao 
  stopifnot(!any(site.chao.suitable & !site.abundance))
  chao <- rep(NA, length(site.abundance))
  chao[site.chao.suitable] <- 
    tapply(diversity$Measurement[diversity$Is_Chao_suitable], 
           droplevels(diversity$SSS[diversity$Is_Chao_suitable]), 
           function(m) sum(m>0) + (((sum(m==1) * (sum(m==1)-1)) / (2*(sum(m==2)+1)))))
  return (chao)
}

.SiteRarefiedRichness<-function(diversity, site.abundance){
  .Log("Computing Rarefied Species Richness\n")
  # Compute rarefied richness only for studies where Diversity_metric_is_suitable_for_Chao 
  # is TRUE and all measurements within the study are integers.
  
  # http://stackoverflow.com/questions/3476782/how-to-check-if-the-number-is-integer
  # TODO all.equal
  study.rsr.suitable <- 
    tapply(diversity$Measurement[diversity$Diversity_metric_is_suitable_for_Chao],
           diversity$SS[diversity$Diversity_metric_is_suitable_for_Chao], 
           function(m) all(floor(m)==m))
  
  # Studies for which Diversity_metric_is_suitable_for_Chao is FALSE will have 
  # a value of NA in study.rsr.suitable
  study.rsr.suitable[is.na(study.rsr.suitable)] <- FALSE
  site.rsr.suitable <- study.rsr.suitable[tapply(diversity$SS, diversity$SSS, unique)]
  
  diversity$Is_RSR_suitable <- site.rsr.suitable[diversity$SSS]
  
  # Only abundance metrics should be marked as suitable for rarefied richness 
  stopifnot(!any(site.rsr.suitable & !site.abundance))
  
  rsrich <- rep(NA, length(site.abundance))
  
  if (any(diversity$Is_RSR_suitable)){
    siteTotalAbund<-aggregate(Measurement~SSS,data=droplevels(
      diversity[diversity$Is_RSR_suitable,]),FUN=sum)
    siteTotalAbund$SS<-diversity$SS[match(siteTotalAbund$SSS,diversity$SSS)]
    studyMinSiteAbund<-aggregate(Measurement~SS,data=siteTotalAbund,FUN=
                                   function(x) min(x[x>0]))
    diversity$studyMinSiteAbund<-studyMinSiteAbund$Measurement[match(
      diversity$SS,studyMinSiteAbund$SS)]
    
    sample.n <- split(diversity$studyMinSiteAbund[diversity$Is_RSR_suitable], 
                      droplevels(diversity$SSS[diversity$Is_RSR_suitable]))
    values <- split(diversity$Measurement[diversity$Is_RSR_suitable],
                    droplevels(diversity$SSS[diversity$Is_RSR_suitable]))
    
    rsr<-matrix(nrow=length(which(site.rsr.suitable)),ncol=1000)
    for (i in 1:1000){
      rsr[,i]<-mapply(function(vals,samp){
        sp<-as.character(1:length(vals))
        ind<-rep(sp,vals)
        if (length(ind)==0){
          sr<-0
        } else {
          rare.sp<-sample(ind,samp[1])
          sr<-length(table(rare.sp))
        }
        return(sr)
      },values,sample.n)
    }
    
    rsrich[site.rsr.suitable]<-apply(rsr,1,mean)
    
  }
  
  
  
  
  return(rsrich)
  
}
