##### Compare different data.
require(pacman)
p_load(tidyverse)
source("data.R")

# set variable names.
data.var = 'data'
features.var = 'features'
stim.var = 'Stimulation_treatment'
site.var = "Image_Metadata_Site"
track.var = 'track_id'
Grouping = "Grouping"
exp.var = "exp"

# load pre saved data.
load("data/data.three.RData")
load("data/features.three.RData")

source("load_multiple.R")

# function call first source function below.
#experiment_comparison_plots(data.three, features.three, stim.var, site.var, track.var, exp.var)

#' Compares trends per feture between treatment to negative control, excluding the positive control
#' It returns two plots, one a zscore for each siRNA per experiment plate and one averaged for for each plate.
#'
#' @param data.three data frame containg the raw time series data. With experiment identifier.
#' @param features.three  data frame containing the calculated features. With experiment identifier.
#' @param stim.var column name for the treatment identifiers.
#' @param site.var name of the site identifier, e.g. track_id
#' @param track.var name of the track/trajectorie identifier, e.g. track_id
#' @param exp.var name of the experiment identifier
#'
#' @return
#' @export
#'
#' @examples
experiment_comparison_plots <- function(data.three, features.three, stim.var, site.var, track.var, exp.var) {
  

  comb_p = lapply(data.three[[exp.var]] %>% unique, function(k){
    ###  separate per plate--------- 
    dt.data = data.three %>% filter(get(exp.var) == k) %>% ungroup
    features = features.three %>% filter(get(exp.var) == k) %>% ungroup
    
    # separate into -CTRL and rest without +CTRL. Pool different wild types.
    f.ex.neg = features %>% filter(!get(stim.var) %like% "CTRL") %>% mutate(Stimulation_treatment = case_when(get(stim.var) %like% "WT" ~ "WT", TRUE ~ Stimulation_treatment))
    f.neg = features %>% filter(get(stim.var) %like% "\\-CTRL") 
    # calc mean and sd of -CTRL per plate
    z_m = f.neg  %>% summarise_if(is.numeric,mean,na.rm = T) %>% select(-site.var, -track.var, -exp)
    z_sd = f.neg  %>% summarise_if(is.numeric,sd,na.rm = T) %>% select(-site.var, -track.var, -exp)
    names = f.ex.neg %>% pull(get(stim.var))
    # z score formula.
    all = (f.ex.neg %>% select(-stim.var, -track.var, -site.var, -exp.var) - (z_m %>% slice( rep(1,nrow(f.ex.neg))))) / (z_sd %>% slice( rep(1,nrow(f.ex.neg)))) 
    all[[stim.var]] = names
    z_per_si = all %>% group_by(.dots=c(stim.var)) %>% summarise_if(is.numeric, mean, na.rm = T)
    z_per_si$plate = paste("Plate ", k)
    z_per_si
  })
  # bind everyting together and convert to long.
  a = do.call(rbind,comb_p)
  a = a %>% gather(key = "feature", value = zscore, - Stimulation_treatment,-plate)
  # discard features we dont use anymore.
  a  = a %>% filter(!feature %in% c("baseline.after.slope"))
  
  ### plotting ------
  plot1 = ggplot(a, aes(y = feature, x = forcats::fct_rev(as.factor(get(stim.var))), fill = zscore)) + 
    facet_grid(plate~., scale = "free_y") + 
    geom_raster() + ggplotTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    #geom_text(aes(label = round(zscore,2))) +  
    labs(y = "Feature", x = "siRNA Treatment", fill = "z-score")
  ggsave(plot1, filename= "plotsforthesis/final/zscorepertreatment.svg", device = "svg", width = 15.8, height = 9.19)
  
  ### summarised plot ------
  plot2 = ggplot(a %>% group_by(plate, feature) %>% summarise(zscore = mean(zscore, na.rm = T)), aes(x = feature, y = forcats::fct_rev(as.factor(plate)), fill = zscore)) + 
    #facet_grid(plate~., scale = "free_y")+ 
    geom_raster() + ggplotTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    #geom_text(aes(label = round(zscore,2))) + 
    labs(y = "Experiment", x = "Feature", fill = "z-score") + labs_pubr(base_size = 20)
  plot(plot1)
  plot(plot2)
  ggsave(plot2, filename= "plotsforthesis/final/zscoreperplate.svg", device = "svg", width = 15.8, height = 9.19)
  return(a %>% group_by(plate, feature) %>% summarise(zscore = mean(zscore, na.rm = T)))
}

experiment_comparison_per_ctrl(data.three, features.three, exp.var, stim.var, Grouping, site.var, track.var)

#' Z score per CTRL then averaged per feature
#' not used in the thesis
#'
#' @param data.three raw time series data of multiple experiments
#' @param features.three feature data of of multiple experiments
#' @param exp.var  column identifier of experiment identifier.
#' @param stim.var column identifier stimulation_treatment
#' @param Grouping Grouping per CTRL
#' @param site.var column identifier of site
#' @param track.var column identifier of the track variable
#'
#' @return None
#' @export
#'
#' @examples
experiment_comparison_per_ctrl <- function(data.three, features.three, exp.var, stim.var, Grouping, site.var, track.var) {

  comb_p = lapply(data.three[[exp.var]] %>% unique, function(k){
    ## case of one plate. ---
    dt.data = data.three %>% filter(get(exp.var) == k) %>% ungroup
    features = features.three %>% filter(get(exp.var) == k) %>% ungroup
    
    unique(dt.data %>% filter(get(stim.var) %like% "CTRL") %>% pull(stim.var))
    
    #separate CTRL ad
    #dt.ctrl = dt.data %>% filter(get(stim.var) %like% "CTRL")
    #f.ex.ctrl = features %>% filter(!get(stim.var) %like% "CTRL")
    
    nr_groups = length(unique(dt.data %>% pull(Grouping)))
    comb = lapply(seq(nr_groups) , FUN = function(i){
      group = paste0("\\-CTRL ", (i -1))
      z_m = features %>% filter(get(stim.var) %like% group) %>% summarise_if(is.numeric,mean,na.rm = T) %>% select(-site.var, -track.var)
      # frame all but 
      z_sigma = features %>% filter(get(stim.var) %like% group) %>% summarise_if(is.numeric,sd,na.rm = T) %>% select(-site.var, -track.var)
      m_per_si =  features %>%
        filter(!get(stim.var) %like% group) %>% filter(!get(stim.var) %like% "CTRL") %>%
        group_by(.dots = stim.var) %>% summarise_if(is.numeric,mean,na.rm = T) %>% select(-site.var, -track.var, stim.var)
      
      names = m_per_si %>% pull(stim.var)
      m_per_si = m_per_si %>% select(-stim.var, -exp.var)
      z_m = z_m %>% select(-exp.var)
      z_sigma = z_sigma %>% select(-exp.var)
      z_per_si = ((m_per_si - (z_m %>% slice( rep(1,nrow(m_per_si)))) ) / (z_sigma %>% slice( rep(1,nrow(m_per_si)))))
      #z_per_si[[stim.var]] = names
      r = z_per_si %>% summarise_all(mean, na.rm = T)
      r = r %>% mutate_if(is.infinite, function(y){NaN})
    })
    a = do.call(rbind,comb)
    a = a %>% mutate(CTRL = paste0("CTRL ", seq(nr_groups) -1), plate = k) 
    a
  })
  comb_p = do.call(rbind,comb_p)
  
  comb_p = comb_p %>% gather(key = "feature", value = zscore, - CTRL, -plate)
  ggplot(comb_p, aes(x = feature, y = forcats::fct_rev(as.factor(CTRL)), fill = zscore)) + 
    facet_grid(plate~., scale = "free_y")+ 
    geom_raster() + ggplotTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    #geom_text(aes(label = round(zscore,2))) + 
    labs(y = "Experiment", x = "Feature", fill = "z-score") + labs_pubr(base_size = 20)
  
  ggplot(comb_p %>% group_by_(.dots = c("CTRL", "plate", "feature")) %>% summarise(zscore = mean(zscore)), aes(x = feature, y = forcats::fct_rev(as.factor(plate)), fill = zscore)) + 
    #facet_grid(plate~., scale = "free_y")+ 
    geom_raster() + ggplotTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    #geom_text(aes(label = round(zscore,2))) + 
    labs(y = "Experiment", x = "Feature", fill = "z-score") + labs_pubr(base_size = 20)
}
