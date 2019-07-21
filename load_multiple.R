#'  @Author Cédric Walker


# load some data
if(FALSE){
source("data.R")
source("normalize_data_frame.R")

experiments = c("Coralie/NIH3T3/siPOOLs/20190408_systII_siPOOLs_plate1_and_2_singlePulse_changed_order/20190408_150050_815",
                "Coralie/NIH3T3/siPOOLs/20190311_systII_siPOOLs_plate1_and_2_singlePulse/20190311_102308_731",
                "Coralie/NIH3T3/siPOOLs/20190401_systII_siPOOLs_plate1_and_2_singlePulse/20190401_155035_686")
data_comb = load_multiple(experiments = experiments, mnt = "/run/user/1000/gvfs/afp-volume:host=izbhelsinki,volume=imaging.data")
data_combined = data_comb$data
features.combinded = data_comb$features
features.n.combined = data$features.n
save(data_combined, "data/data_combinded")
save(features.combinded, "data/features.combinded")
save(features.n.combined, "data/features.n.combined")
exp.var = "exp"
site.var = "Image_Metadata_Site"
track.var = "track_id"
time.var = "Image_Metadata_T"
stim.var = "Stimulation_treatment"
erk.ratio.var = "erk.ratio" 
source("onedclustering.R")
features.excluding.pos.ctrl = features.n.combined %>% filter(!stim.var %like% "\\+CTRL")
cluster_1D(features.excluding.pos.ctrl, feature = "amplitude.max" , stim.var = stim.var, track.var = track.var, exp.var = exp.var)
cluster = cluster_1D(features.clean.pooled, feature = "amplitude.max" , stim.var = stim.var, track.var = track.var, exp.var = exp.var)
plot_1D_clustering(cluster$freq.table, palette = "Set1", feature = "amplitude.max")
averages = representatives_1D_clustering(cluster.id = cluster$cluster.id, data = data.three, track.var = track.var, site.var = site.var, time.var = time.var, erk.ratio.var = erk.ratio.var, exp.var = exp.var)
plot_representatives_1D_clustering(representative.table = averages$lines_per_cl, feature = "amplitude.max", palette = "Set1", )
}

#' Load multiple experiments at once. Safes features and data to ./data/
#'
#' @param experiments list of experiment names
#' @param mnt mounting point. is different between operating systems. Single Letter with colon on Windows, Ubuntu is default.
#' Mac I dont know but smthing with /Volume/
#' Filter information is saved to ./data/
#' 
#' @param ... keyword variable pairs for load_data function
#'
#' @return combined tables containing all experiments in a list. , data, features, features.normalized.
#' @export
#'
#' @examples
load_multiple = function(experiments, mnt,...){
  source("data.R")
  source("normalize_data_frame.R")
  data1 = load_data(mnt = mnt,
                   experiment = experiments[1], ...)
  i = 1
  # safe experiment desc as exeriment variable, touple(exp.var, site.var, track.var) is unique per trajectoriy
  # data$data[[exp.var]] = regmatches(experiment,
  #                                   regexpr("(?<=/)\\d+.+(?=/)", experiment, perl = T))
  # data$features[[exp.var]] = regmatches(experiment,
  #                                       regexpr("(?<=/)\\d+.+(?=/)", experiment, perl = T))
  data = data1$data
  features = data1$features
  
  features.n = normalize_frame(data1$features, -site.var, -track.var, -exp.var)
  
  for(experiment in experiments[-1]){
    i = i +1 
    datatmp = load_data(mnt = mnt,
                         experiment = experiment,
                         ...)
    data = rbind(data, datatmp$data)
    features = rbind(features, datatmp$features)
    features.n = rbind(normalize_frame(datatmp$features, -site.var, -track.var, -exp.var), features.n)
    
  }
  # TODO filter information
    return(list(data = data, features = features, features.normalized = features.n))
}

#' Pooles different controls and WT together
#'
#' @param features feature table
#' @param stim.var variable string of column name with siRNA name
#'
#' @return Table, no features are changed!
#' @export
#'
#' @examples
pool_treatments = function(features, stim.var){

  f.p = features %>% mutate(st = case_when(
    get(stim.var) %like% "\\+CTRL" ~ "+CTRL",
    get(stim.var) %like% "\\-CTRL" ~ "-CTRL",
    get(stim.var) %like% "WT" ~ "WT",
    # there are some experiements with multiple siRNA Wells, 
    # Stimulation_treatment %like% "PP2A" ~ "PP2A", 
    # Stimulation_treatment %like% "DUSP4" ~ "DUSP4",
    # Stimulation_treatment %like% "DUSP6" ~ "DUSP6",
    # Stimulation_treatment %like% "DUSP26" ~ "DUSP26",
    # Stimulation_treatment %like% "DUSP22" ~ "DUSP22",
    # Stimulation_treatment %like% "DUSP9" ~ "DUSP9",
    # Stimulation_treatment %like% "DUSP4" ~ "DUSP4",
    # Stimulation_treatment %like% "RKIP" ~ "RKIP",
    TRUE ~ as.character(get(stim.var))
  ) )  %>% select(-stim.var) %>%  rename_(.dots = setNames("st", stim.var))
  return(f.p)
}
