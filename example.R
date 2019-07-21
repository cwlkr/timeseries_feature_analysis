# exmaple workflow

source("data.R")
source("normalize_data_frame.R")
source("onedclustering.R")     
source("load_multiple.R")

# path to experiment on the izbhelsinki
experiments = c("Coralie/NIH3T3/siPOOLs/20190408_systII_siPOOLs_plate1_and_2_singlePulse_changed_order/20190408_150050_815",
                "Coralie/NIH3T3/siPOOLs/20190311_systII_siPOOLs_plate1_and_2_singlePulse/20190311_102308_731",
                "Coralie/NIH3T3/siPOOLs/20190401_systII_siPOOLs_plate1_and_2_singlePulse/20190401_155035_686")
# only run once per exp setup. saves whole compination of experiments
exp.mult = load_multiple(experiments)
exp.mult.data = exp.mult$data
exp.mult.features.n = exp.mult$features.n # features normalized, per experiment
exp.mult.features = exp.mult$features
save(exp.mult.data, "data.RData")
save(exp.mult.features.n, "features.n.RData")
save(exp.mult.features, "features.RData")

# name of loaded element is the same as saved object
load("data.RData")
load("features.n.RData")
load("features.RData")
#set column information
exp.var = "exp"
site.var = "Image_Metadata_Site"
track.var = "track_id"
time.var = "Image_Metadata_T"
stim.var = "Stimulation_treatment"
erk.ratio.var = "erk.ratio" 


# plate comparison plot
source("zscoretoctrl.R")

experiment_comparison_plots(exp.mult.data, exp.mult.features, stim.var, site.var, track.var, exp.var) 
  

# pool +/-ctrls and wt and remove +CTRL 
features.n.pooled = pool_treatments(exp.mult.features.n, stim.var) %>% filter(!get(stim.var) %like% "\\+CTRL")

source("feature_correlation.R")
plot_feature_correlation(features.n.pooled, track.var, site.var, exp.var, T)

#set the feature for the clustering
feature = "amplitude.max"

# load group information
functional_f = readxl::read_xlsx(path = "grouped_clusterings/siPOOLs_grouping.xlsx", sheet = 1)
f_na = functional_f[,6] %>% slice(4:7) %>% t()  %>% unname %>% c
f_nr = functional_f[,5] %>% slice(4:7) %>% t()  %>% unname %>% as.numeric() %>% c
functional_f %>% slice(3) %>% select(1:2) -> names
functional_f %>% slice(4:n()) %>% select(1:2) -> functional
colnames(functional) = names
group.names = data.frame(Groups = f_nr, name = f_na)


# sets for each trajectory the corresponding cluster
clustering = cluster_1D(features.n.pooled , feature, stim.var, track.var, exp.var, Grouping = functional)
# plotting the cluster proportions for each siRNA treatment.
plot_1D_clustering(clustering$freq.table, facet = T, group.names = group.names, feature = feature)
# plot cluster average trajectories
rep.n =representatives_1D_clustering(clustering$cluster.id, exp.mult.data, track.var, site.var, time.var, erk.ratio.var, exp.var, Grouping = functional)
plot_representatives_1D_clustering(rep.n, feature = feature, group.names = group.names)
