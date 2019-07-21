#' @author Cédric Walker 

#' Add the baseline features to the TS plot
#'
#' @param baseline output of the baseline function
#' @param series_len lenghts of the TS
#'
#' @return None
#' @export
#'
#' @examples
add_baseline_to_plot = function(baseline, series_len){
  lines(1:10, rep(baseline$before, 10), col = "red")
  lines((series_len-9):series_len, rep(baseline$after, 10), col = "red")
}
#' Adds the decay to the plot
#'
#' @param series raw time series 
#' @param lambda decay coefficient
#' @param max.amp output of peak detection
#'
#' @return None
#' @export
#'
#' @examples
add_dec_to_plot = function(series, lambda, max.amp){
  lines(c(max.amp$time.max,length(series)), series[max.amp$time.max] + lambda * c(0, length(series) - max.amp$time.max), col = "yellow")
  
}
#' adds the growth to the plot
#'
#' @param series raw time series 
#' @param lambda growth coefficient
#' @param max.amp output of peak detection
#'
#' @return None
#' @export
#'
#' @examples
add_growth_to_plot = function(series, lambda, max.amp){
  x = c(0,50)
  reg = (( x )*  lambda)  - ((lambda * max.amp$time.max) - series[max.amp$time.max])
  lines(x, reg, col = "yellow")
}

#' Add dip to plot
#'
#' @param dip_amp dip object
#' @param series raw time series.
#'
#' @return None
#' @export
#'
#' @examples
add_dip_to_plot = function(dip_amp, series){
  x = rep(dip_amp$time.min,10)
  y1 = series[dip_amp$time.min]
  y =sample(c(y1, y1 + dip_amp$min), 10, replace = T)
  lines(x,y, col = "green")
}
#' add amplitude to plot
#'
#' @param max_amp output of peak detection
#' @param series raw time series.
#'
#' @return None
#' @export
#'
#' @examples
add_amp_to_plot = function(max_amp, series){
  x = rep(max_amp$time.max,10)
  y1 = series[max_amp$time.max]
  y =sample(c(y1, y1 - max_amp$max), 10, replace = T)
  lines(x,y, col = "green")
}
#' Add the FWHM to the plot
#'
#' @param fmwh output of FWHM wrapper function of the single time series detection
#' @param series raw time series 
#'
#' @return None 
#' @export
#'
#' @examples
fwhm_to_plot = function(fwhm, series){
  
  x = c(fwhm$left, fwhm$right)
  dec = fwhm$right - as.integer(fwhm$right)
  y = rep( (series[as.integer(fwhm$right)] * (1 -dec)) + (  series[as.integer(fwhm$right) +1 ] * dec) ,2 )
  lines(x, y, col = "blue")
  
}

#' Title
#' Connot be looped as r plots use lazy paradigma
#'
#' @param slope.after.max output of slope after maximum function
#' @param max.amp output of peak detection function
#' @param series raw time series
#'
#' @return None
#' @export
#'
#' @examples
plot_slopes_after = function(slope.after.max, max.amp, series){
  slope.after.max.n = length(slope.after.max)
  ampm = max.amp$time.max
  c = round(length(series[max.amp$time.max:length(series)])/slope.after.max.n)
  i = 0
  b = series[ampm+((i*c) +(c/2))] - (slope.after.max[[i+1]] * (ampm+((i*c) +(c/2))))
  x = c(ampm + (i*c), ampm + c + (i*c))
  y = slope.after.max[[i+1]] * x + b
  lines(x,y, col = "red")
  i = 1
  b = series[ampm+((i*c) +(c/2))] - (slope.after.max[[i+1]] * (ampm+((i*c) +(c/2))))
  x = c(ampm + (i*c), ampm + c + (i*c))
  y = slope.after.max[[i+1]] * x + b
  lines(x,y, col = "red")
  i = 2
  b = series[ampm+((i*c) +(c/2))] - (slope.after.max[[i+1]] * (ampm+((i*c) +(c/2))))
  x = c(ampm + (i*c), ampm + c + (i*c))
  y = slope.after.max[[i+1]] * x + b
  lines(x,y, col = "red")
  # i = 3
  # b = series[ampm+((i*c) +(c/2))] - (slope.after.max[[i+1]] * (ampm+((i*c) +(c/2))))
  # x = c(ampm + (i*c), ampm + c + (i*c))
  # y = slope.after.max[[i+1]] * x + b
  # lines(x,y, col = "red")
}


#' plot a time series with its features
#'
#' @param series raw time series
#' @param baseline output of baseline function
#' @param dip output of dip detection function
#' @param amp output of amp detection function
#' @param fwhm output of fwhm extraction function
#' @param dec output of decay extraction function
#' @param growth output of growth extraction function
#' @param slope.after.max output of function
#'
#' @return None
#' @export
#'
#' @examples
series_plot_feat = function(series, baseline, dip, amp, fwhm, dec, growth, slope.after.max, stim.times){
  plot(series, type = "l", xlab = 'Time (min)', ylab = 'ERK-KTR Cytoplasmic to Nuclear Intensity', ylim=c(0,1.2))
  add_baseline_to_plot(baseline, series_len = length(series))
  add_dip_to_plot(dip, series)
  add_amp_to_plot(amp, series)
  fwhm_to_plot(fwhm, series)
  add_dec_to_plot(series, lambda = dec,max.amp = amp)
  add_growth_to_plot(series, lambda = growth,max.amp = amp)
  plot_slopes_after(slope.after.max, amp, series)
  rug(stim.times, col = "red", lw = 2)
}

# -- Plotting stuff
if(FALSE){
# dt.data %>% group_by(.dots = c("Image_Metadata_Site",stim.var, "objNuc_TrackObjects_Label")) %>%
#   mutate(nuc_erk = get(cyto.erk)/get(nuc.erk)) %>%
#   do(cbind(., FeatAllFeat(.$nuc_erk, basal = 0.5, 1, 20))) %>% View(.)

# pdf("PP2A.pdf")
# 
# data = dt.data %>% unite("unique" , c(objNuc_TrackObjects_Label, Image_Metadata_Site), sep = "_") %>% filter(Stimulation_treatment == "PP2A")
# 
# 
# for(i in unique(data$unique)){
#   single_timeseries_coralie = data %>%  filter(unique == i) %>%
#     mutate(nuc_erk = objCyto_ring_Intensity_MeanIntensity_imErk/objNuc_Intensity_MeanIntensity_imErk) %>% select(nuc_erk) %>% melt() %>% pull()
# 

  baseline = series_baseline(single_timeseries_coralie, 10)
  max_amp = FeatMaxAmplitude(single_timeseries_coralie, basal = baseline$before)
  fmwh = FeatFWHM(single_timeseries_coralie,  basal = baseline$before)
  dip_amp = dip_amplitude(single_timeseries_coralie, stim.time = 10, baseline)
  dec = FeatHalfMaxDec(single_timeseries_coralie)
  growth = FeatHalfMaxGrow(single_timeseries_coralie)
  sla = afterpeakdecay(single_timeseries_coralie, max_amp)
  series_plot_feat(single_timeseries_coralie, baseline, dip_amp, max_amp, stim.time, fmwh, dec, growth, sla)
  single_timeseries_coralie_treat = data %>% filter(unique == k[i]) %>%
    select(Stimulation_treatment) %>% pull()

  title(as.character(single_timeseries_coralie_treat)[1])

# 
# }
# dev.off()
# 
}
# ### end plot stuff ------