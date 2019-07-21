#'  @Author Cédric Walker


# justification of randomForest oulier detection

#' Keeps track of filtering steps in the data loading process. And creates tables reporting proportaions for all filters,
#' and for the random Forest tree filtering separate. saving to ./data/...filename
#' function is called from the load_data function of data.R
#'
#' @param prediction vector for each trajectory indicating 0 for positive and 1 for negative
#' @param dt.data.w.features_class features with NA replaced with 1000
#' @param pos.ctrl.features features for trajectories of the +CTRL
#' @param classifier the random forest classifier
#' @param rec_filter output of receptor intensity filtering
#' @param max_filter output of maximum filtering
#' @param total_bef_filter the total number of trajectories before filtering.
#' @param stim.var column identifier of the stimulation_treatment.
#' @param filename where to write filter data information.
#'
#' @return list of tables containg filter information of a data frame. 
#' @export
#'
#' @examples
validate_outliers = function(prediction,
                             features.na.replaced,
                             pos.ctrl.features,
                             classifier,
                             rec_filter,
                             max_filter,
                             total_bef_filter,
                             stim.var,
                             filename){
  
  #----------
  require(tidyverse)
  require(randomForest)
  require(data.table)

  outliers = which(prediction == 1)
  a = table(features.na.replaced[outliers,]$Stimulation_treatment)
  b = table(features.na.replaced[-outliers,]$Stimulation_treatment) # remove lots of plus control which is expected!
  
  
  pos.ctrl.features_class = pos.ctrl.features
  pos.ctrl.features_class[is.na(pos.ctrl.features_class)] = 1000
  prediction_pos = predict(classifier, pos.ctrl.features_class)
  outliers_pos = which(prediction_pos == 1)
  c = table(pos.ctrl.features_class[outliers_pos,]$Stimulation_treatment)
  d = table(pos.ctrl.features_class[-outliers_pos,]$Stimulation_treatment)
  sum(table(pos.ctrl.features_class[outliers_pos,]$Stimulation_treatment))/length(pos.ctrl.features_class$Stimulation_treatment)
  
  at = data.frame(a) %>% rename("Freq1" = Freq)
  bt = data.frame(b) %>% rename("Freq2" = Freq)
  
  outtable = dplyr::full_join(at, bt) %>% mutate(Freq1 = case_when(is.na(Freq1) ~ 0, TRUE ~ as.double(Freq1)))
  outtable_tree = outtable %>%
    mutate(total = (Freq2+Freq1), perc_outlier = Freq1/(Freq2+Freq1)) %>% mutate(perc_outlier = if_else(!is.na(perc_outlier),round(perc_outlier, digits = 3),0))
  outtable_tree = `colnames<-`(outtable_tree, c("Treatment", "outlier", "not outlier", "total", "precentage of outliers"))
  
  outtable_gen = full_join(outtable %>% rename("Treatment" = Var1), rec_filter %>% rename_(.dots = setNames(stim.var, "Treatment"))) %>%
    full_join(.,max_filter %>% rename_(.dots = setNames(stim.var, "Treatment"))) %>% select(-Freq2) %>% 
    full_join(total_bef_filter  %>% rename_(.dots = setNames(stim.var, "Treatment")))
  
  outtable_gen[is.na(outtable_gen)] = 0
  outtable_gen = outtable_gen %>% mutate(perc_filtered = (Freq1 + rec_filter + max_filter)/total, traj_left = total - (Freq1 + rec_filter + max_filter) ) %>%
    mutate(perc_filtered = round(perc_filtered, digits = 3)) %>% rename("tree_filter" = Freq1)
  
  return( list(
    tree_validation = outtable_tree %>% arrange(desc(get('precentage of outliers'))),
    general_validation = outtable_gen %>% arrange(desc(perc_filtered))
  ))
  #-----------
}




#' Plot a subset of 50 preidcted as positive and 50 predicted as negative for visual control.
#'
#' @param all.ex.pos.features.clean fetures without pos CTRL predicted as pos
#' @param features.outliers features predicted as negative. 
#' @param dt.data  raw timeseries data, unfiltered. 
#' Can not be the return of load_data, needs to be the manuell step until random forest filtering.
#'
#' @return None
#' @export
#'
#' @examples
plot_test = function(all.ex.pos.features.clean, features.outliers, dt.data){
  pdf("bs_forest/validate_outliers_clean.pdf")
  all.ex.pos.clean.unique = all.ex.pos.features.clean %>% unite(unique, c(Image_Metadata_Site, track_id, exp)) %>% select(unique)
  set.seed(156)
  k = sample(all.ex.pos.clean.unique$unique,50)
  clean.data.to.plot = dt.data %>% unite(unique, c(Image_Metadata_Site, track_id, exp)) %>% filter(unique %in% k)
  for(i in unique(clean.data.to.plot$unique)){
    single_timeseries_coralie = clean.data.to.plot %>%  filter(unique == i) %>%
      mutate(nuc_erk = objCyto_ring_Intensity_MeanIntensity_imErk/objNuc_Intensity_MeanIntensity_imErk) %>% select(nuc_erk) %>% melt() %>% pull()
    
    baseline = series_baseline(single_timeseries_coralie, 10)
    max_amp = FeatMaxAmplitude(single_timeseries_coralie, basal = baseline$before)
    fmwh = FeatFWHM(single_timeseries_coralie,  basal = baseline$before)
    dip_amp = dip_amplitude(single_timeseries_coralie, stim.time = 10, baseline)
    dec = FeatHalfMaxDec(single_timeseries_coralie)
    growth = FeatHalfMaxGrow(single_timeseries_coralie)
    slope.after.max = afterpeakdecay(single_timeseries_coralie, max_amp, 4)
    series_plot_feat(single_timeseries_coralie, baseline, dip_amp, max_amp, stim.time, fmwh, dec, growth, slope.after.max,4)
    single_timeseries_coralie_treat = clean.data.to.plot %>% filter(unique == k[i]) %>%
      select(Stimulation_treatment) %>% pull()
    
      title(as.character(single_timeseries_coralie_treat)[1])
    
  }
  dev.off()
  #outliers
  features.outliers
  pdf("bs_forest/validate_outliers_outliers.pdf")
  
  all.ex.pos.outliers.unique = features.outliers %>% unite(unique, c(Image_Metadata_Site, track_id)) %>% select(unique)
  set.seed(156)
  k = sample(all.ex.pos.outliers.unique$unique,50)
  outlier.data.to.plot = dt.data %>% unite(unique, c(Image_Metadata_Site, track_id)) %>% filter(unique %in% k)
  for(i in unique(outlier.data.to.plot$unique)){
    single_timeseries_coralie = outlier.data.to.plot %>%  filter(unique == i) %>%
      mutate(nuc_erk = objCyto_ring_Intensity_MeanIntensity_imErk/objNuc_Intensity_MeanIntensity_imErk) %>% select(nuc_erk) %>% melt() %>% pull()
    
    baseline = series_baseline(single_timeseries_coralie, 10)
    max_amp = FeatMaxAmplitude(single_timeseries_coralie, basal = baseline$before)
    fmwh = FeatFWHM(single_timeseries_coralie,  basal = baseline$before)
    dip_amp = dip_amplitude(single_timeseries_coralie, stim.time = 10, baseline)
    dec = FeatHalfMaxDec(single_timeseries_coralie)
    growth = FeatHalfMaxGrow(single_timeseries_coralie)
    slope.after.max = afterpeakdecay(single_timeseries_coralie, max_amp, 3)
    series_plot_feat(single_timeseries_coralie, baseline, dip_amp, max_amp, stim.time, fmwh, dec, growth, slope.after.max)
    single_timeseries_coralie_treat = data %>% filter(unique == k[i]) %>%
      select(Stimulation_treatment) %>% pull()
    
    title(as.character(single_timeseries_coralie_treat)[1])
    
  }
  dev.off()
  
}

