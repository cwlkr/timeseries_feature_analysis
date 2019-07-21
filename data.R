#'  @Author Cédric Walker

if(!require(pacman)){
  install.packages("pacman")
}

#rm(list  = ls())
pacman::p_load(data.table, ggplot2, tidyverse, readxl, randomForest, TSexploreR, R.utils)


#' All purpose function for loading a single- or multipulse timelapse experiment.
#' Filters accoring to max erk, receptor intensity and a random forest algorithm for detecting failures in feature calculations.
#' Saves data, features and filtering information to subolderer /data in currend wd (default - can be changed for data and features), . 
#'
#'
#' @param nuc.erk column identifier of nuclear erk mean intensity
#' @param cyto.erk column identifier of cytosolic erk mean intensity
#' @param time.var column identifier of the time varable
#' @param stim.var column identifier of the Stimulation treatment variable
#' @param stim.time.var  column identifier for the stimulation time variable in the metadata file
#' @param group.var vector of column names which is a unique touple per trajectory.
#' @param erk.ratio.var "new variable the erk ratio will be saved in. Default is  erk.ratio
#' @param meta.grouping column identifier in the metadata file for the grouping of mean plots
#' @param track.var column identifier for tracking identifier
#' @param site.var column identifier for site identifier 
#' @param data.folder folder name the data is saved to, is added to the infered experiment name. full path or relative path to working directory, without /, e.g., c:/data
#' @param experiment folder the experiment data is in. (Default loads experiment 20190408)Needs to point to root foolder where folder output-data and receptor are in.
#' @param plot.filename (default infers name from experiment name)if plot = T plots will be saved accoring to this filename path.
#' @param plot boolean variable. Deafult is FALSE. Create mean-plots per Grouping.
#' @param write.all boolean variable. Default is TRUE. Saves data and features as csv files to data subfolder in current cwd.
#' @param erk_a (deafult is reasonable)quantile threshold for maxumum erk filtering. 
#' @param rec_a (deafult is reasonable)quantile threshold for recepter intensity level filtering. Lower Threshold. 
#' @param rec_b (deafult is reasonable)quantile threshold for recepter intensity level filtering. Upper Threshold.
#' @param normalize boolean variable. Default is FALSE. Set TRUE if tracjectories should be normalized before feature extraction. Normalize threshold accoring to the average of the first 5 time points.
#'
#' @return a list with the filtered data set and feature calculations as well as documentation of evefy filtering step.
#' @export
#'
#' @examples load_data()
    load_data = function(
    nuc.erk = "objNuc_Intensity_MeanIntensity_imErk",
    cyto.erk = "objCyto_ring_Intensity_MeanIntensity_imErk",
    #cyto.erk = "objCyto_Intensity_MeanIntensity_imErk" # alterntive in older Experiments
    time.var = "Image_Metadata_T",
    stim.var = "Stimulation_treatment",
    stim.time.var = "Stimulation_time",
    group.var = c("Image_Metadata_Site", "track_id"),
    erk.ratio.var = "erk.ratio",
    meta.grouping = "Grouping",
    track.var = "track_id",
    site.var = "Image_Metadata_Site",
    exp.var = "exp",
    
    data.folder = NULL,
    
    
    #laod dataset
    ## example single ------
    mnt =  "/run/user/1000/gvfs/afp-volume:host=izbhelsinki,volume=imaging.data", # ubuntu
    # mnt = 'y:'                                                               # windows
    # experiment = "Coralie/NIH3T3/siPOOLs/20181210_systIII_siPOOLs_plate1_multipulses_II/20181210_075541_575/Merged"
    # experiment = "Coralie/NIH3T3/siPOOLs/20190221_systIII_siPOOLs_plate1_singlePulses/20190221_090712_551/Merged"
    experiment = "Coralie/NIH3T3/siPOOLs/20190408_systII_siPOOLs_plate1_and_2_singlePulse_changed_order/20190408_150050_815",
    plot.filename = NULL,
    plot = F,
    write.all = T,
    erk_a = 0.99,
    rec_a = 0.05,
    rec_b = 0.95,
    normalize = F
){
  ## set path information ----------
  metadatapath = paste(mnt, experiment, sep = "/")
  path = paste(mnt, experiment, "output-data", sep = "/") # or just the path to tCoursesSelected.csv/ Merged Folder
  
  ## infere filename information based on experiment name if filename has not been specified manually--------
  name = regmatches(experiment,
                    regexpr("(?<=/)\\d+.+(?=/)", experiment, perl = T))

  # Read the csv files of the raw data and receptor data ---------
  s.file = list.files(path = paste(mnt, experiment, "output-data", sep = "/"), pattern = "tCoursesSelected.csv.gz", recursive = TRUE, full.names = TRUE)
  print(s.file)
  dt.data = fread(file = s.file)
  r.file = list.files(path = paste(mnt, experiment, "output-data", sep = "/"), pattern = "(rec)(.)*.csv.gz", recursive = TRUE, full.names = TRUE)
  print(r.file)
  
  receptordata = fread(r.file)
  
  #load metadata and process metadata information ----------
  m.file = list.files(path = metadatapath, pattern =  "plate(.*).xlsx" , recursive = FALSE, full.names = TRUE) 
  print(m.file)
  metadata = as.data.table(read_xlsx(m.file))[-(1:2)]
  metanamecol = paste0(as.vector(metadata[1,]), "") # dirty hack to get "real" vector
  metadata = `colnames<-`(metadata, metanamecol)[-1,]
  
  # calculate and add erk ratio -----------
  dt.data = dt.data[, erk.ratio := get(cyto.erk)/get(nuc.erk)]
  
  # extract stimulation time(s) from the metadata
  stim.times <- metadata %>% select(contains(stim.time.var), stim.var) %>% 
    select(contains(stim.time.var)) %>% slice(1)  %>% c(., recursive=TRUE) %>% as.numeric()
  
  
  #just to make sure that grouping and stim_treatment are ok! ---------
  groups <- metadata %>% # group_by(.dots = meta.grouping) %>% 
    select(stim.var, meta.grouping, Position, Well) %>% rename("Image_Metadata_Site" = "Position", "Metadata_Well" = "Well") %>%
    mutate("Grouping" = as.numeric(get(meta.grouping)), Image_Metadata_Site = as.numeric(Image_Metadata_Site), Metadata_Well = as.numeric(Metadata_Well)) 
  
  dt.data =   left_join(dt.data, groups)
  
  # set up filter gathering information. For each treatment, count number of trajectories
  total_bef_filter = dt.data %>% group_by(.dots = c(track.var, site.var)) %>% summarise(st = unique(get(stim.var))) %>%
    rename_(.dots=setNames("st" ,stim.var)) %>% ungroup %>% group_by(.dots = stim.var) %>% summarise(total = n())
  
  ### max erk ratio filtering-------
  # Check for abnormal max.
  max_m_erk = dt.data %>% group_by(.dots = c(group.var)) %>% summarise(max_erk = max(get(erk.ratio.var), na.rm = T)) %>% pull(max_erk) %>% max()
  # plot histogram for manuall inspection
  dt.data %>% group_by(.dots = c(group.var)) %>% summarise(mean_erk = max(get(erk.ratio.var), na.rm = T)) %>% pull(mean_erk) %>% hist(500, freq = F)
  abline(v = dt.data %>% group_by(.dots = c(group.var)) %>% 
           summarise(max_erk = max(get(erk.ratio.var), na.rm = T)) %>% pull(max_erk) %>% quantile(erk_a) , col = "red")
  # abline(v = dt.data %>% group_by(.dots = c(group.var)) %>% 
  #          summarise(max_erk = mean(get(erk.ratio.var), na.rm = T)) %>% pull(max_erk) %>% mean() )

  x = dt.data %>% group_by(.dots = c(group.var)) %>% summarise(max_erk = max(get(erk.ratio.var), na.rm = T)) %>% filter(max_erk < quantile(max_erk,erk_a)) %>% select(-max_erk)
  dt.mean.reduced = left_join(x, dt.data)
  
  # process and save max_filter information to keep track of number of trajectories filtered.
  max_filter = anti_join(dt.data %>% select(track.var, site.var, stim.var), dt.mean.reduced %>% select(track.var, site.var, stim.var)) %>% 
    group_by(.dots = c(track.var, site.var)) %>% summarise(st = unique(get(stim.var))) %>%
    rename_(.dots=setNames("st" ,stim.var)) %>% ungroup %>% group_by(.dots = stim.var) %>% summarise(max_filter = n())
  max_filter_id = anti_join(dt.data %>% select(track.var, site.var, stim.var), dt.mean.reduced %>% select(track.var, site.var, stim.var)) %>% 
    group_by(.dots = c(track.var, site.var)) %>% summarise(st = unique(get(stim.var))) %>%
    rename_(.dots=setNames("st" ,stim.var)) %>% ungroup %>% select(-stim.var)
  dt.data = dt.mean.reduced
  
  ### recpetor filtering -----------
  # plot receptor intensity histogram
  cut = (receptordata$obj_Rec_Intensity_MeanIntensity_imRecCorrOrig)
  breaks =  seq(0,1000,1) / (1000 / max(receptordata$obj_Rec_Intensity_MeanIntensity_imRecCorrOrig))
  p = hist(cut, freq = F, breaks = breaks, main = "Receptor Intensity per TS", xlab = "Receptor Intensity")
  
  abline(v=quantile(receptordata$obj_Rec_Intensity_MeanIntensity_imRecCorrOrig,rec_a, na.rm = T), col = "red")
  abline(v = quantile( receptordata$obj_Rec_Intensity_MeanIntensity_imRecCorrOrig,rec_b, na.rm = T) , col = "red")
  
  ci_l = quantile( receptordata$obj_Rec_Intensity_MeanIntensity_imRecCorrOrig,rec_a, na.rm = T)
  ci_u = quantile( receptordata$obj_Rec_Intensity_MeanIntensity_imRecCorrOrig,rec_b, na.rm = T)
  receptordata %>% filter(obj_Rec_Intensity_MeanIntensity_imRecCorrOrig > ci_l &
                            obj_Rec_Intensity_MeanIntensity_imRecCorrOrig < ci_u) %>% select(site.var, track_id) -> reduce
  dt.rec.reduced = inner_join(reduce, dt.data)

  # keep track of trajectories filtered by the receptor intensity filtering.
  rec_filter = anti_join(dt.data %>% select(track.var, site.var, stim.var), reduce %>% select(track.var, site.var)) %>% 
    group_by(.dots = c(track.var, site.var)) %>% summarise(st = unique(get(stim.var))) %>%
    rename_(.dots=setNames("st" ,stim.var)) %>% ungroup %>% group_by(.dots = stim.var) %>% summarise(rec_filter = n())
  rec_filter_id = anti_join(dt.data %>% select(track.var, site.var, stim.var), reduce %>% select(track.var, site.var)) %>% 
    group_by(.dots = c(track.var, site.var)) %>% summarise(st = unique(get(stim.var))) %>%
    rename_(.dots=setNames("st" ,stim.var)) %>% ungroup %>% select(-stim.var)
  
  dt.data = dt.rec.reduced
  
  # normalize data on the first 5 data points. Feature extraction is then made on normalized trajectories.
  if(normalize){
    dt.data.m = dt.data %>% filter(get(time.var) < 5) %>% group_by(.dots = c(site.var, track.var)) %>% summarise(diff = 1 - mean(get(erk.ratio.var)))
    dt.data = left_join(dt.data, dt.data.m) %>% mutate(erk.ratio.n = get(erk.ratio.var) + diff) %>% select(-erk.ratio.var) %>%
      rename_(.dots = setNames("erk.ratio.n", erk.ratio.var))
   #print(paste("normalizeing", as.character(dt.data %>% pull(get(erk.ratio.var)) %>% summary())) )
  }
  ####  plot average trajectories accoring to groups---------------
  
  if(plot){
    source('~/MasterProjects/timeseries_analysis/plots.R')
    if(is.null(plot.filename)){
      plot.filename = paste0("data/mean_plot", regmatches(experiment,
                                                          regexpr("(?<=/)\\d+.+(?=/)", experiment, perl = T)))
    }
    plot_average_per_group(dt.data = dt.data,plot.filename = pdf.filename, stim.times = stim.times, meta.grouping = meta.grouping,
                           time.var = time.var, erk.ratio.var = erk.ratio.var, stim.var = stim.var,
                           vlines = F)
  }
  
  # calc features and predict useless features
  if(length(stim.times) == 1){
  source('timeseries_analysis_functions_single.R')
  dt.data.w.features <- run_feature_extraction_single(data = dt.data, stim.time = stim.times, group.var, erk.ratio.var, track.var = track.var, site.var = site.var, stim.var = stim.var) 
  
  classifier = readRDS("bs_forest/bs_forest.classifier_data20190311_plate1and2")
  
  all.ex.pos.features = dt.data.w.features %>% filter(!(Stimulation_treatment %like% "\\+CTRL"))
  pos.ctrl.features = dt.data.w.features %>% filter((Stimulation_treatment %like% "\\+CTRL"))
  
  
  dt.data.w.features_class = all.ex.pos.features
  dt.data.w.features_class[is.na(dt.data.w.features_class)] = 1000
  prediction = predict(classifier, dt.data.w.features_class)
  length(prediction[which(prediction == 1)])/length(prediction)
  
  ## check what stim_treatments get discarded!
  outliers = which(prediction == 1)
  
  
                                  
  all.ex.pos.features.clean = all.ex.pos.features[-outliers,]
  features.outliers = all.ex.pos.features[outliers,]
  features.clean = rbind(all.ex.pos.features.clean, pos.ctrl.features)
  features.clean[[exp.var]] = name
  dt.data[[exp.var]] = name
  } # else{
  #   source("timeseries_analysis_functions_multi.R")
  #   
  # }
  # generate filter information tables and validation sets of the random forest calculations.
  source("outlierdetection_validation.R")
  validation = validate_outliers(prediction,
                    dt.data.w.features_class,
                    pos.ctrl.features,
                    classifier, rec_filter,
                    max_filter,
                    total_bef_filter,
                    stim.var,
                    regmatches(experiment,
                               regexpr("(?<=/)\\d+.+(?=/)", experiment, perl = T)))
  
  if(write.all){
    if(is.null(data.folder)){
      filename.features = paste0("data/","features_", name, ".csv")
      filename.data = paste0("data/", "data_", name, ".csv")
      filename.xlsx = paste0("data/", "filtering_table_",name,".xlsx")
    }else{
      filename.features = paste0(data.folder, "/features_",name, ".csv")
      filename.data = paste0(data.folder, "/data_",name, ".csv")
      filename.xlsx = paste0(data.folder, "filtering_table_",name,".xlsx")
    }
    
    require(openxlsx)
    list_of_datasets <- list("tree_validation" = validation$tree_validation, "general_validation" = validation$general_validation)
    write.xlsx(list_of_datasets, file =filename.xlsx)
    fwrite(features.clean, filename.features)
    fwrite(dt.data, filename.data)
  }
  return(list(data = dt.data, features = features.clean, validation = validation, filter_data = list(max_filter = max_filter_id, rec_filter = rec_filter_id, tree_filter = features.outliers %>% select(site.var, track.var))))
}

## Is not excecuted. For loading default vales of column names into memory.    
if(FALSE){
  nuc.erk = "objNuc_Intensity_MeanIntensity_imErk"
  cyto.erk = "objCyto_ring_Intensity_MeanIntensity_imErk"
  #cyto.erk = "objCyto_Intensity_MeanIntensity_imErk" # alterntive in older Experiments
  time.var = "Image_Metadata_T"
  stim.var = "Stimulation_treatment"
  stim.time.var = "Stimulation_time"
  group.var = c("Image_Metadata_Site", "track_id")
  erk.ratio.var = "erk.ratio"
  meta.grouping = "Grouping"
  track.var = "track_id"
  site.var = "Image_Metadata_Site"
  
  filename.data = "tdata_20190408_systIII_siPOOLs_plate1and2_singlePulses"
  
  
  filename.features = "features_20190408_systIII_siPOOLs_plate1and2_singlePulses"
  
  #laod dataset
  ## example single ------
  mnt =  "/run/user/1000/gvfs/afp-volume:host=izbhelsinki,volume=imaging.data" # ubuntu
  # mnt = 'y:'                                                               # windows
  # experiment = "Coralie/NIH3T3/siPOOLs/20181210_systIII_siPOOLs_plate1_multipulses_II/20181210_075541_575/Merged"
  # experiment = "Coralie/NIH3T3/siPOOLs/20190221_systIII_siPOOLs_plate1_singlePulses/20190221_090712_551/Merged"
  experiment = "Coralie/NIH3T3/siPOOLs/20190408_systII_siPOOLs_plate1_and_2_singlePulse_changed_order/20190408_150050_815"
  pdf.filename = "plots_multi_20190311.pdf"
  plot = F
  write.all = T
  erk_a = 0.99
  rec_a = 0.05
  rec_b = 0.95
}
    
#' helper function for creating a ggplot theme wiht pubr titles
#'
#' @param in.font.base 
#' @param in.font.axis.text 
#' @param in.font.axis.title 
#' @param in.font.strip 
#' @param in.font.legend 
#'
#' @return ggplot theme
#' @export
#'
#' @examples
ggplotTheme = function(in.font.base = 12,
                       in.font.axis.text = 12,
                       in.font.axis.title = 12,
                       in.font.strip = 14,
                       in.font.legend = 12) {
  library(ggpubr)
  loc.theme =
    theme_bw(base_size = in.font.base, base_family = "Helvetica") +
    theme(
      panel.spacing = unit(1, "lines"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.25),
      axis.text = element_text(size = in.font.axis.text),
      axis.title = element_text(size = in.font.axis.title),
      strip.text = element_text(size = in.font.strip, face = "bold"),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = in.font.legend),
      legend.key.height = unit(1, "lines"),
      legend.key.width = unit(2, "lines")) + labs_pubr() + theme_pubclean(base_size = 14)
}
  