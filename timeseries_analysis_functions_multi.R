pacman::p_load(TSexploreR, data.table, tidyverse)



#' calculate the first order partial derivative of the timeseries 
#'
#' @param series raw time series vector
#'
#' @return time series vector first order partial derivatives
#' @export
#'
#' @examples
calc_derivatives = function(series){
  return(series - lag(series))
}

#' calculate the TS baseline before and after
#'
#' @param series raw time series data
#' @param stim.time numerical, fist occurance of stimulation
#' @param robust boolean, if True use median, else (default) mean
#'
#' @return named list with baseline before and after as well as standart deviation
#' @export
#'
#' @examples
series_baseline_mp = function(series, stim.time, robust = F){
  before = mean(series[1:stim.time])
  after = mean(series[(length(series) - 9) : (length(series))])
  if(robust){
    before = median(series[1:stim.time])
    after = median(series[(length(series) - 9) : (length(series))])
  }
  
  
  return(list(before = before, 
              after  = after,
              diff = abs(after - before),
              before.sd = sd(series[1:stim.time]),
              after.sd = sd(series[(length(series) - 9) : (length(series))])
  ))
}

#' calculate and locate dip of TS
#'
#'
#' @param series raw time series as vector
#' @param stim.time time of first stimulation
#' @param baseline output of series_baseline_mp function
#'
#' @return named list with dip amplitude and time
#' @export
#'
#' @examples
dip_amplitude_mp <- function(series, stim.time, baseline){
  # Shift the data to have basal at 0, this is important so that maximum peak is
  # measured in amplitude instead of absolute value
  series <- series - baseline$before
  if(  min(series[stim.time:which.max(series)]) > 0){dip = 0
  } else{dip = min(series[(stim.time + 1):which.max(series)])} # the real dip cannot already happen at stim time either
  return(list(min = abs(dip), time.min = stim.time -1 + which.min(series[stim.time:which.max(series)])))
}

#' response delay, as diff between each peak and the simulations
#'
#' @param stim.times stimulation times
#' @param max.amp output of mp_amp_wrapper
#'
#' @return vector, delay for each stimulation to peak
#' @export
#'
#' @examples
response_delay = function(stim.times, max.amp){
  return( list(delay = abs(as.numeric(names(max.amp$extrema)) - stim.times)))
}


#' wrepper of TSexploreR::MPFeatTrend_extrema, with some assumptions about the TS
#' 
#' @param series raw time series vector
#' @param stim.time time of first stimulation
#' @param baseline output of series_baseline_mp function
#' @param robust boolean, if True (default) median regreesion, else mean regression
#' @param basewin basewindow, area of series that are considered baseline before and after the stimulation are disregared during feature extraction
#'
#' @return named list of peaks and their time, model used as well as the amplitude difference between the peaks
#' @export
#'
#' @examples
mp_amp_wrapper = function(series, stim.time, baseline, robust = T, basewin = 10){
  
  red = series[basewin : (length(series) - basewin)]
  
  amp = TSexploreR::MPFeatTrend_extrema(red, 5, "maxi", robust)
  extrema = amp$extremes
  x.pos  = as.numeric(names(extrema)) + basewin - 1
  names(extrema) <- x.pos
  amplitude = extrema - baseline$before
  if (robust){
    l.m = mblm(extrema~x.pos)
  }else {
    l.m = lm(extrema~x.pos)
  }
  if(!is_empty(which(amplitude < 0))){
    extrema = extrema[-which(amplitude < 0)]
    amplitude = amplitude[-which(amplitude < 0)]
  }
  return( list(trend = l.m$coefficients[2],
               model = l.m,
               extrema = extrema,
               amplitude = amplitude,
               diff = abs((extrema - baseline$before) - lag(extrema - baseline$before))[-1]))
}

#' load all featues to data frame
#'
#' @param series raw time series vector
#' @param stim.times stimulation timms 
#' @param basewin basewindow, area of series that are considered baseline before and after the stimulation are disregared during feature extraction
#'
#' @return data.frame of all features
#' @export
#'
#' @examples
all_features_df_multi = function(series, stim.times, basewin){
  
  bline = series_baseline_mp(series, stim.times[1])
  dip_amplitude = dip_amplitude_mp(series, stim.times[1], bline)
  am = mp_amp_wrapper(series, baseline = bline, basewin = basewin)
  delay = response_delay(stim.time = stim.times, am)
  
  return(data.frame(c("baseline" = bline,
               "amplitude" = am$amplitude,
               "amplitude.times" = am$extrema,
               "amplitude.coeff" = am$trend,
               amplitude.diff = am$diff,
              "dip" = dip_amplitude,
              "delay" = delay)))
  }


run_feature_extraction_multi = function(data, stim.times, group.var, erk.ratio.var, track.var, site.var, stim.var){
  dt.w.features <- data %>%
    group_by_at(.vars = group.var)  %>% 
    do(all_features_df_multi(.[[erk.ratio.var]], stim.times, 10)) 
  stim.df = data %>% select(.dots = c(site.var, track.var, stim.var)) %>% distinct()
  colnames(stim.df) = c(site.var, track.var, stim.var)
  return(left_join(dt.w.features, stim.df))
}

#' Looping through lines does not work because of lazy r base plot evaluation
#' plotting it as single connecting line. Does not work for each multipeak TS currently
#'
#' @param max.amp output of mp_amp_wrapper
#'
#' @return None
#' @export
#'
#' @examples
add_amp_to_plot = function(max.amp){
  x = rep(as.numeric(names(max.amp$amplitude[1])) , 2)
  y = c(max.amp$extrema[1], max.amp$extrema[1] - max.amp$amplitude[1])
  for(a in  2:(length(max.amp)) ){
    x = c(x, rep(as.numeric(names(max.amp$amplitude[a])) , 2))
    y = c(y, c(max.amp$extrema[a], max.amp$extrema[a] - max.amp$amplitude[a]))
  }
  lines(x, y, col = "green")
}

if(FALSE){
  dt.data.mp = fread(file = "tdata_20181206_systIII_siPOOLs_plate1_multipulses_I.csv")
  series =  dt.data.mp %>% filter(objNuc_TrackObjects_Label == 15 & Image_Metadata_Site == 15) %>%
    mutate(nuc_erk = objCyto_ring_Intensity_MeanIntensity_imErk/objNuc_Intensity_MeanIntensity_imErk) %>% select(nuc_erk) %>% melt() %>% pull()
  
bline = series_baseline_mp(series, stim.times[0])
plot(series, type = "l", xlab = "Time (min)", ylab = "ERK-KTR C/N")
add_baseline_to_plot(bline, series_len = length(series)) # located in plot_features_per_series.R
add_dip_to_plot(dip_amplitude(series, stim.time = 10, bline), series) # located in plot_features_per_series.R
rug(stim.times, col = "red", lw = 2)
basewin = stim.times[1]
am = mp_amp_wrapper(series, baseline = bline, basewin = basewin, robust = T)
abline(am$model$coefficients, col = "blue")
add_amp_to_plot(am)
# add amplitude to multipeak features. Not loopable beacause of r base plot lazy eval
##
x = rep(as.numeric(names(am$amplitude[1])) , 2)
y = c(am$extrema[1], am$extrema[1] - am$amplitude[1])
lines(x, y, col = "green")
x = rep(as.numeric(names(am$amplitude[2])) , 2)
y = c(am$extrema[2], am$extrema[2] - am$amplitude[2])
lines(x, y, col = "green")
x = rep(as.numeric(names(am$amplitude[3])) , 2)
y = c(am$extrema[3], am$extrema[3] - am$amplitude[3])
lines(x, y, col = "green")
x = rep(as.numeric(names(am$amplitude[4])) , 2)
y = c(am$extrema[4], am$extrema[4] - am$amplitude[4])
lines(x, y, col = "green")
#add_amp_to_plot(am)
group.var = c("Image_Metadata_Site", "objNuc_TrackObjects_Label")

run_feature_extraction_multi(dt.data.mp, stim.times = stim.times, group.var,erk.ratio.var = "erk.ratio", track.var =  "objNuc_TrackObjects_Label", stim.var = "Stimulation_treatment",site.var = site.var)


}