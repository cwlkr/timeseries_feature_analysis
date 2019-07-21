# for singlepulse and non noisy signal.
require(tidyverse)
require(TSexploreR)
require(data.table)

#' Calcu
#'
#' @param series as a numerical vector
#'
#' @return fist discrete derivate of the time series
#' @export
#'
#' @examples
calc_derivatives = function(series){
  return(series - lag(series))
}
 
#' calculate baseline features of a single pulse time series.
#'
#' @param series as a vector 
#' @param stim.time time the stimulation occurs
#' @param basewindow window the baseline is consiered. Default is 10 frames
#' @param robust boolean. If T median is used, default is mean.
#'
#' @return list of baseline features. Mean and sd
#' @export
#'
#' @examples
series_baseline = function(series, stim.time, basewindow = 10, robust = F){
  basewindow = basewindow -1
  before = mean(series[1:stim.time])
  after = mean(series[(length(series) - basewindow) : (length(series))])
  if(robust){
    before = median(series[1:stim.time])
    after = median(series[(length(series) - basewindow) : (length(series))])
  }
  after_y = series[(length(series) - basewindow) : (length(series))]
  fit = lm(y~after.slope, data.frame(after.slope = 1:length(after_y), y = after_y))
  return(list(before = before, 
              after  = after,
              diff = (before - after),
              before.sd = sd(series[1:stim.time]),
              after.sd = sd(series[(length(series) - 9) : (length(series))]),
              after.slope = coef(fit)[2]
         ))
}

#' calculates fwhm features as a wrapper of the walking algorithm of the TSexplorer package.
#'
#' @param series raw time series as a vector.
#' @param basal baseline before information
#' @param dip.amp amplitude of the dip, calculated by the dip_amplitude function
#' @param basewindow window of the baseline, reduces noise by using baseline as a constraint.
#'
#' @return list of FWHM features.
#' @export
#'
#' @examples
fwhm_wrapper = function(series,  basal, dip.amp, basewindow = 10){
  red = series[dip.amp$time.min :(length(series) - basewindow -1)]
  fwhm = FeatFWHM(y = red, basal = basal)
  
  return( list(
    fwhm = fwhm$fwhm,
    right = fwhm$right + dip.amp$time.min -1,
    left = fwhm$left + dip.amp$time.min -1
  ))
}
# only possible if stim.time is known and max amplitude measured correctly

#' Title
#'
#' @param series time series as a vector
#' @param stim.time stimulation time
#' @param baseline baseline features
#'
#' @return amplitude of the dip und time of the amplitude dip.
#' @export
#'
#' @examples
dip_amplitude <- function(series, stim.time, baseline){
  # Shift the data to have basal at 0, this is important so that maximum peak is
  # measured in amplitude instead of absolute value
  series <- series - baseline$before
  red_y = series[(stim.time +1):   (which.max(series[stim.time: (length(series) - stim.time)])  +stim.time - 1 ) ]  
  offset = stim.time
  if(  min(red_y) > 0){
    dip = 0
  } else{
    dip = min(red_y)
    } # the real dip cannot already happen at stim time either
  return(list(min = abs(dip), time.min =offset + which.min(red_y)))
}

#' fitting an exponential decay to the slope. Does not work well.
#'
#' @param series time series as a vector
#' @param stim.time stimulation time 
#'
#' @return exponential decay parameter lamda
#' @export
#'
#' @examples
exp_decay = function(series, stim.time){
  series = single_timeseries_coralie
  peak = FeatMaxAmplitude(series)
  max.time = peak$time.max
  #follow until d/dt shifts
  y = series[-(1:max.time)]
  min = min(y)
  max = max(y)
  min.time = which.min(y) + max.time - 1
  return (log(max/min)/abs(min.time - max.time))
}

#' Delay of the response calculated as tje dofferemce between time of amplitude max and stimulation time.
#'
#' @param stim.time stime the stimulation occurse
#' @param max.amp  return value of the max_amp_wrapper function. or named vector with field time.max.
#'
#' @return named vector with the delay value
#' @export
#'
#' @examples
response_delay = function(stim.time, max.amp){
  return( c("delay" = max.amp$time.max -stim.time)) 
}
#' Wrapper of the TSexploreR::FeatMaxAmplitude function, uses the constrained that there can be no peak before the simtulation occured.
#'
#' @param series vector containing time series values.
#' @param stim.time numerical. Number of frame the stimulation occured.
#' @param baseline output of the series_baseline function.
#' @param dip.amp  output of the dip_amp function.
#'
#' @return named vector containing the amplitude value and the frame it was detected.
#' @export
#'
#' @examples
max_amp_wrapper = function(series, stim.time, baseline, dip.amp){
  max_amp = FeatMaxAmplitude(series[dip.amp$time.min:(length(series) - stim.time)], basal = baseline$before)
  if (max_amp$max < 0) { max = NA}else{max = max_amp$max}
  return( list(
    max = max,
    time.max = (max_amp$time.max + dip.amp$time.min - 1) 
  ))
}

#' Calculates the growth as a line between the left point of the FWHM and the amplitude max
#'
#' @param series vector containing time series values.
#' @param max.amp output of the max_amp_wrapper function.
#' @param fwhm output of the fwhm_wrapper function.
#'
#' @return named vector conaining the slope of the growh or na if fwhm was not found.
#' @export
#'
#' @examples
growth_half_max= function(series, max.amp, fwhm){
  if(!is.na(fwhm$left)){
  f_p = fwhm$left - as.integer(fwhm$left)
  y_fwhm_inter = (series[fwhm$left] * (1-f_p)) + (series[fwhm$left + 1] * (f_p))
  growth = (series[max.amp$time.max] - y_fwhm_inter)/(max.amp$time.max - fwhm$left)
  return(c(growth = growth))
  } else {
    return (NA)
    }
}

#' Calculates the decay as a line between the right point of the FWHM and the amplitude max
#'
#' @param series vector containing time series values.
#' @param max.amp output of the max_amp_wrapper function.
#' @param fwhm output of the fwhm_wrapper function.
#'
#' @return named vector conaining the slope of the decay or na if fwhm was not found.
#' @export
#'
#' @examples
decay_half_max = function(series, max.amp, fwhm){
  if(!is.na(fwhm$right)){
    f_p = fwhm$right - as.integer(fwhm$right)
    y_fwhm_inter = (series[fwhm$right] * (1-f_p)) + (series[fwhm$right + 1] * (f_p))
    decay = (y_fwhm_inter - series[max.amp$time.max])/(fwhm$right - max.amp$time.max )
    return(c(decay = decay))
  } else {
    return (NA)
  }
}

#' For time series slope calculation. Splits the time series into n section an calculated the slope as a linear regression.
#'
#' @param series vector containing time series values.
#' @param max.amp output of the max_amp_wrapper function.
#' @param n numerical number of sections
#'
#' @return named vector conaining each the slope coefficient for each section.
#' @export
#'
#' @examples
afterpeakdecay = function(series, max.amp, n=3){
  cut = series[max.amp$time.max : length(series)]
  cut_N  = round(length(cut)/n)
  list = lapply(0:(n-1), function(i){
    if(i != (n-1)){
      j = (i*cut_N) + 1
      fit = lm( y~x, data.frame(y =cut[ j:((i+1) * cut_N) ], x =  (j:((i+1) * cut_N) ))) 
      coef(fit)[2]
      }else{
        j = (i*cut_N) + 1
        fit = lm( y~x, data.frame(y =cut[ j:length(cut) ], x =  (j:length(cut) ))) 
        coef(fit)[2]
  }
  })
  list
}


#' Describes the noise of a time series as the change in sign of the first derivative.
#'
#' @param series a vector containing raw time series data
#' @param stim.time number of frame the stimulation occurs.
#'
#' @return numerical. noise measurement. 
#' @export
#'
#' @examples
change_in_derivatives = function(series, stim.time){
  X = calc_derivatives(series)[(stim.time):length(series)][-1]
  last_sign = 1
  sign_changes = 0
 
  for(x in X){
    if (x == 0){
      sign = -1
    } else{
        sign = x / abs(x)
    }
    if(sign == -last_sign){
        # calculate mean and sd in window around the change, norm series in window and scale by diff 
          sign_changes = sign_changes + 1
          last_sign = sign
    }
    last_x = x
  }
  return(sign_changes)
}

#' Title
#'
#' @param series 
#' @param stim.time 
#' 
#' @return
#' @export
#'
#' @examples
all_features_df = function(series, stim.time){

  baseline = series_baseline(series, stim.time)
  dip_amp = dip_amplitude(series, stim.time = stim.time, baseline)
  max_amp = max_amp_wrapper(series, stim.time, baseline, dip.amp = dip_amp)
  fwhm = fwhm_wrapper(series,  basal = baseline$before, dip_amp)
  delay = response_delay(10, max_amp)
  # rates!
  slope.after.max = afterpeakdecay(series, max_amp,3)
  dec = decay_half_max(series, max_amp, fwhm)
  growth = growth_half_max(series, max_amp, fwhm)
  noise = change_in_derivatives(series, stim.time)
  return( data.frame(c("baseline" = baseline,
            'amplitude' = max_amp,
            'FWHM' = fwhm,
            'dip' = dip_amp,
            'delay' = unname(delay),
            'growth' = unname(growth),
            'decay' = unname(dec),
            mean = mean(series),
            sd = sd(series),
            noise = noise,
            slope.after.max = slope.after.max
            )))
}

#' Title
#'
#' @param data time series data. Data frame conaining muliple trajectories in long format.
#' @param stim.time stimulation time, numerical
#' @param group.var vector containing column names wich are sufficent for unique identification of a time series.
#' @param erk.ratio.var column variable in which the time sereis values are encoded.
#' @param track.var column name of the track_id identifier. 
#' @param site.var  column name of the site identifier
#' @param stim.var name the siRNA treatment is encoded in.
#'
#' @return data frame with all extracted features.
#' @export
#'
#' @examples
run_feature_extraction_single = function(data, stim.time, group.var, erk.ratio.var, track.var, site.var, stim.var){
  dt.w.features <- data %>%
    group_by_at(.vars = group.var)  %>% 
    do(all_features_df(.[[erk.ratio.var]], stim.time)) 
  stim.df = data %>% select(.dots = c(site.var, track.var, stim.var)) %>% distinct()
  colnames(stim.df) = c(site.var, track.var, stim.var)
  return(left_join(dt.w.features, stim.df ))
}



