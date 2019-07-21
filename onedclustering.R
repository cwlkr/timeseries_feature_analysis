#'  @Author Cédric Walker


#' Creates a one dimensional clustering based on one z-normalised feature
#'
#' @param feature feature column name
#' @param stim.var column identifier of Stimulation_treatment
#' @param track.var column identifier of the track identifier.
#' @param exp.var column identifier of the experiement identifier.
#' @param Grouping grouping variable, e.g, based on localisation.
#' @param normalized.features.ex.pos features without pos.ctrl normalized for each experiment idividually.
#'
#' @return the frequency table of the clustering as well as the unique id to cluster pair.
#' @export
#'
#' @examples
cluster_1D = function(normalized.features.ex.pos, feature, stim.var, track.var, exp.var = NULL,Grouping = NULL ){
  data.red = normalized.features.ex.pos %>% select(stim.var, track.var, feature, site.var, exp)
  tot = data.red %>%  group_by(.dots = c(stim.var)) %>% summarise(total = n())
  
  f.dt = data.red %>% mutate(cl = case_when(
    get(feature) > 1 ~ "upper",
    get(feature) < -1 ~ "lower",
    TRUE ~ "middle"
  ))
  f.dt.count = f.dt %>% group_by(.dots = c("cl", stim.var)) %>% summarise(count = n()) %>%
    complete(Stimulation_treatment, cl = c("lower", "middle", "upper") , fill = list(count = 0))
  
  x = full_join(f.dt.count, tot) %>% mutate(freq = count /total)
  #fill up empty stuff
  
  if(!is.null(Grouping)){
    x = Grouping %>% rename_(.dots = setNames("Target",stim.var)) %>% right_join(., x)
  }
  return(list( freq.table = x, cluster.id = f.dt))
}

#' Plot the frequencies of the one dimensional clustering
#'
#' @param data freq.table output of cluster_1D. 
#' @param palette color palette for ggplot.
#' @param group.names names of grouping used. e.g. for localisttion c(0="CTRL", 1 = "Adapters")
#' @param facet boolean, if T faceting accoring to Group
#' @param feature name of the clustered feature.
#'
#' @return ggPlot, cluster frequencies plot.
#' @export
#'
#' @examples
plot_1D_clustering = function(freq.table, palette, group.names = NULL, facet = F, feature){
  if(!is.null(freq.table$Group) && !is.null(group.names)){
    group.names_v<- setNames(as.character(group.names$name), as.character(group.names$Groups))
    freq.table = freq.table %>% mutate(Group = if_else(is.na(Group), "Control", dplyr::recode(Group, !!!group.names_v)))
    
  }
  
  gg = ggplot(freq.table,aes(x = Stimulation_treatment, y = count, fill = cl)) + 
    geom_bar(position = position_fill(reverse = TRUE),stat = "identity") +
    scale_y_continuous(labels = scales::percent_format()) +
    
    ggplotTheme() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20)) +
    scale_fill_brewer(palette = palette) +  theme(text = element_text(size=25)) +
    labs(y = "precentage of trajectories", fill = "Cluster", x = "siRNA treatment", title = feature) 
  if(facet){
    if(is.null(freq.table$Group)){
      print("No Group variable found!")
    }else{
      gg = gg + facet_grid(~Group, scale = "free_x", space = "free") 
    }
  }
  plot(gg)
  return(gg)
}

#' creates mean plots per clustering. 
#'
#' @param cluster.id cluster.id ouput of cluster_1D function
#' @param data raw time series data.
#' @param track.var column identifier of the track_id
#' @param site.var column identifier of the site.
#' @param time.var column identifier of the time variable
#' @param erk.ratio.var column identifier of the erk measurement
#' @param exp.var column identifier of the experiment identifier
#' @param Grouping if a grouping variable is set, cluster averages are also calcualted for each group seperately.
#'
#' @return list of two tables containing the average trajectories per cluster, and if possible per Group, else NULL
#' @export
#'
#' @examples
representatives_1D_clustering = function(cluster.id, data, track.var, site.var, time.var, erk.ratio.var, exp.var = NULL, Grouping = NULL){
  
  cluster.id %>% left_join(.,data, by=c(track.var, site.var, exp.var)) %>% group_by(.dots=c("cl", time.var)) %>%
    summarise(erk = mean(get(erk.ratio.var)), n = n()) -> lines_per_cl
  
  lines_per_g_per_cl = NULL
  if(!is.null(Grouping)){
    
    Grouping %>% rename_(.dots = setNames("Target",stim.var)) %>% right_join(., cluster.id) %>% left_join(.,data, by=c(track.var, site.var, exp.var)) %>% 
      group_by(.dots=c("Group", "cl", time.var)) %>% summarise(erk = mean(get(erk.ratio.var)), n = n())  -> lines_per_g_per_cl
  }
  list(lines_per_cl = lines_per_cl,lines_per_g_per_cl= lines_per_g_per_cl )
}


#' Plots the ouput of the representatives_1D_clustering function.
#'
#' @param data ouput of the representatives_1D_clustering function.
#' @param palette color pallete default is Set1
#' @param feature string identifier of the clustered feature.
#' @param group.names names of grouping used. e.g. for localisttion c(0="CTRL", 1 = "Adapters")
#'
#' @return list of ggplot plots. first with per group and group names.
#' @export 
#'
#' @examples
plot_representatives_1D_clustering = function(representative.table, feature, palette = "Set1", group.names = NULL){
  if(!is.null(representative.table[[2]]$Group) && !is.null(group.names)){
    group.names_v<- setNames(as.character(group.names$name), as.character(group.names$Groups))
    representative.table[[2]] = representative.table[[2]] %>% ungroup %>% mutate(Group = if_else(is.na(Group), "Control", dplyr::recode(Group, !!!group.names_v)))
  }
  gg = ggplot(representative.table[[2]], aes(x = Image_Metadata_T, y  = erk, color = cl))+
    facet_grid(~Group) + 
    scale_color_brewer(palette = palette) + 
    geom_line() + theme_bw() +theme(text = element_text(size=20)) +labs(x = "Real Time" , y = "C/N ERK-KTR ratio", color = "Cluster", title = feature)
  plot(gg)
  ggt = ggplot(representative.table[[1]], aes(x = Image_Metadata_T, y  = erk, color = cl))+
    scale_color_brewer(palette = palette) + 
    geom_line() + theme_bw() + theme(text = element_text(size=20)) + labs(x = "Real Time" , y = "C/N ERK-KTR ratio", color = "Cluster", title = feature) 
  plot(ggt)
  return(list(gg, ggt))
  
}

# not run if sourced----------- 
if(FALSE){
  ## EXAMPLE
  source("data.R")
  source("normalize_data_frame.R")
  source("onedclustering.R")
  source("load_multiple.R")
  load("data/features.three.n.RData")
  load("data/data.three.RData")
  load("data/data.three.n.RData")
  exp.var = "exp"
  site.var = "Image_Metadata_Site"
  track.var = "track_id"
  time.var = "Image_Metadata_T"
  stim.var = "Stimulation_treatment"
  erk.ratio.var = "erk.ratio" 
  features.clean.pooled = pool_treatments(features.three.n, stim.var) %>% filter(!get(stim.var) %like% "\\+CTRL")
  
  functional_f = readxl::read_xlsx(path = "grouped_clusterings/siPOOLs_grouping.xlsx", sheet = 1)
  functional_f %>% slice(3) %>% select(1:2) -> names
  functional_f %>% slice(4:n()) %>% select(1:2) -> functional
  colnames(functional) = names
  group.names = data.frame(Groups = c(0,1,2,3), name = c("Kinases", "Phosphatases", "Adaptors, scaffolds, antagonists", "GTPase and ligants" ))
  group.names_v <- setNames(as.character(group.names$name), as.character(group.names$Groups))
  
  feature = "slope.after.max2"
  st = cluster_1D(features.clean.pooled, feature, stim.var, track.var, exp.var, Grouping = functional)
  pl1 = plot_1D_clustering(st$freq.table, "Set1", facet = T, group.names = group.names, feature = feature)
  rep.n =representatives_1D_clustering(st$cluster.id, data.three, track.var, site.var, time.var, erk.ratio.var, exp.var, Grouping = functional)
  pl23.n = plot_representatives_1D_clustering(rep.n, palette = "Set1", feature = feature, group.names = group.names)
  rep =representatives_1D_clustering(st$cluster.id, data.three.n , track.var, site.var, time.var, erk.ratio.var, exp.var, Grouping = functional)
  pl23 = plot_representatives_1D_clustering(rep, palette = "Set1", feature = feature, group.names = group.names)
  
  ggsave(pl1 + labs_pubr(base_size = 31), filename = paste0("plotsforthesis/final/", feature,"_dist.svg") , device = "svg", width = 40, height = 20 ,units = "in")
  ggsave(pl23.n[[1]], filename = paste0("plotsforthesis/final/", feature,"_mean_per_f_norm.svg") , device = "svg", width = 20, height = 10 ,units = "in")
  ggsave(pl23.n[[2]], filename = paste0("plotsforthesis/final/", feature,"_mean_norm.svg") , device = "svg", width = 12, height = 10 ,units = "in")
  ggsave(pl23[[1]], filename = paste0("plotsforthesis/final/", feature,"_mean_per_f.svg") , device = "svg", width = 20, height = 10 ,units = "in")
  ggsave(pl23[[2]], filename = paste0("plotsforthesis/final/", feature,"_mean.svg") , device = "svg", width = 12, height = 10 ,units = "in")
  
  # --------------
  
  ## for decay and slope after max2 2D ! ---------------------
  # decay and slope per amp cluster?
  feature1 = "amplitude.max"
  feature2 = "slope.after.max2"
  f_a_clust =  cluster_1D(features.clean.pooled, feature1, stim.var, track.var, exp.var, Grouping = functional)
  comb.freq = lapply(c("lower", "middle", "upper"), function(cl_i){
    f_group = f_a_clust$cluster.id %>% filter(cl == cl_i) %>% select(site.var, exp.var, track.var) %>% left_join(.,features.clean.pooled)
    clust = cluster_1D(f_group, feature2, stim.var, track.var, exp.var, Grouping = functional)
    clust$freq.table$cl_1 = cl_i
    clust$freq.table
  })
  comb.cl_id = lapply(c("lower", "middle", "upper"), function(cl_i){
    f_group = f_a_clust$cluster.id %>% filter(cl == cl_i) %>% select(site.var, exp.var, track.var) %>% left_join(features.clean.pooled)
    clust = cluster_1D(f_group, feature2, stim.var, track.var, exp.var, Grouping = functional)
    clust$cluster.id$cl_1 = cl_i
    clust$cluster.id
  })
  do.call(rbind, comb.freq) -> freq
  do.call(rbind, comb.cl_id) -> clust_id
  freq = freq %>% ungroup %>% mutate(Group = if_else(is.na(Group), "Control", dplyr::recode(Group, !!!group.names_v)))
  cln = c(lower = paste0(feature1, ".lower"), middle = paste0(feature1, ".middle"), upper = paste0(feature1, ".upper"))
  freq = freq %>% ungroup %>% mutate(cl_1 =  dplyr::recode(cl_1, !!!cln))
  
  d2_dist = ggplot(freq,aes(x = Stimulation_treatment, y = count, fill = cl)) + 
    geom_bar(position = position_fill(reverse = TRUE),stat = "identity") +
    scale_y_continuous(labels = scales::percent_format()) +
    facet_grid(~cl_1, scale = "free_x", space = "free_x") + 
    ggplotTheme() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
    scale_fill_brewer(palette = "Set1") + 
    labs(y = "precentage of trajectories", fill = "Cluster", x = "siRNA treatment", title = feature2) 
  
  ggsave(d2_dist + labs_pubr(base_size = 14), filename = paste0("plotsforthesis/final/", feature2 , "_per_", feature1, "_dist.svg"), width = 25, height = 10 ,units = "in")
  cl_lines = clust_id %>% left_join(data.three %>% select(track.var, site.var, exp.var, erk.ratio.var, time.var), by = c( "track_id", "Image_Metadata_Site", "exp")) %>% group_by(cl_1, cl, Image_Metadata_T) %>% summarise(erk = mean(get(erk.ratio.var)), n = n())
  cl_lines.n = clust_id %>% left_join(data.three.n %>% select(track.var, site.var, exp.var, erk.ratio.var, time.var), by = c( "track_id", "Image_Metadata_Site", "exp")) %>% group_by(cl_1, cl, Image_Metadata_T) %>% summarise(erk = mean(get(erk.ratio.var)), n = n())
  
  # ---- sort per lower per cond.
  cl_lines = cl_lines %>% ungroup %>% mutate(cl_1 =  dplyr::recode(cl_1, !!!cln))
  cl_lines.n = cl_lines.n %>% ungroup %>% mutate(cl_1 =  dplyr::recode(cl_1, !!!cln))
  
  d2_line = ggplot(cl_lines, aes(x = Image_Metadata_T, y  = erk, color = cl)) +
    facet_grid(~cl_1) +ggplotTheme() + 
    scale_color_brewer(palette = "Set1") + theme_minimal() + 
    geom_line() + labs(x = "Real Time" , y = "C/N ERK-KTR ratio", color = "Cluster", title = feature)
  ggsave(d2_line, filename = paste0("plotsforthesis/final/", feature2, "_per_", feature1, "_line.svg") , device = "svg", width = 20, height = 10 ,units = "in")
  
  d2_lines.n = ggplot(cl_lines.n, aes(x = Image_Metadata_T, y  = erk, color = cl)) +
    facet_grid(~cl_1) +
    scale_color_brewer(palette = "Set1") + theme_minimal() +  
    geom_line() + labs(x = "Real Time" , y = "C/N ERK-KTR ratio", color = "Cluster", title = paste0(feature, " normalized Plot"))
  ggsave(d2_lines.n, filename = paste0("plotsforthesis/final/", feature2, "_per_", feature1, "_line.n.svg") , device = "svg", width = 20, height = 10 ,units = "in")

  
  # investigate individual cluster
  set.seed(616)
  upp_id =st$cluster.id %>% filter(cl = "upper")
  upp_id[sample(1:nrow(upp_id), size = 100),] %>% 
    left_join(data.three , by = c("track_id", "Image_Metadata_Site", "exp")) %>% unite(., "new",  c(track_id, Image_Metadata_Site, exp), sep = "_") -> exp_data
  ggplot(exp_data, aes(x = Image_Metadata_T, y = erk.ratio, group = new)) +
    facet_wrap(~new) +
    geom_line() 
  
}





