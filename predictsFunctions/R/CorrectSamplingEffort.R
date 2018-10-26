CorrectSamplingEffort <-
function(diversity) {
  # Assume that any missing values of sample effort mean equal sampling 
  # effort was applied.
  missing <- is.na(diversity$Sampling_effort)
  cat('Correcting', sum(missing), 'missing sampling effort values\n')
  diversity$Sampling_effort[missing] <- 1
  
  # TODO Check this logic with Tim
  # Rescale sampling effort to have a maximum of 1 in each source
  cat('Rescaling sampling effort\n')
  diversity$Sampling_effort <- do.call('c', tapply(diversity$Sampling_effort, 
                                                   diversity$SS, function(se) return (se/max(se))))
  
  # Correct for sensitivity to sampling effort
  sensitive <- diversity$Diversity_metric_is_effort_sensitive
  cat('Correcting', sum(sensitive), 'values for sensitivity to sampling', 
       'effort\n')
  diversity$Measurement[sensitive] <- diversity[sensitive,'Measurement'] / 
    diversity[sensitive,'Sampling_effort']
  return (diversity)
}
