
#' performes hierarchical clustering on the data set.
#'
#' @param data feature table
#' @param stim.var column name of tretments as string
#' @param site.var site column identifier
#' @param track.var track column identifier
#' @param time.var time column identifier (e.g., Image_Metadata_T)
#' @param k  Number of clusters
#' @param remove.var remove feature 
#' @param filename saves heatmap png to this location
#' @param cluster.col boolean. cluster feature coluumns as well.
#'
#' @return heatmap and dendrogram
#' @export
#'
#' @examples
create_clustering = function(data, stim.var, site.var, track.var , time.var, k, remove.var = NA, filename = NULL, cluster.col = T){
  require(data.table)
  require(tidyverse)
  require(heatmap3)
  data.pooled.n = data
  names = data.pooled.n %>% pull(var = stim.var) 
  cluster.data.n = data.pooled.n %>% ungroup() %>% 
    select(c(-site.var, -track.var, -remove.var)) #, -dip.time.min, -mean, -sd, -FMWH.left, -FMWH.right, -noise) %>% ... = features to deselect!!!
  
  names = cluster.data.n %>% pull(var = stim.var)
    
  cluster.data.n = cluster.data.n %>% select_if(is.numeric) 
  m = as.matrix(cluster.data.n)
  rownames(m) = names
  distance_m = dist(m, method = "euclidean")
  hr <- hclust(distance_m, method = "ward.D2")
  
  ####Create distance matrx and cluster with hc clust
  
  mycl <- cutree(hr, k = k)
  cl = as.factor(mycl)
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  mycolhc <- gg_color_hue(k)
  mycolhc_c <- mycolhc[as.integer(as.vector(cl))] 
  m_clipped = m # hard clipping, ouliers disturbe colors
  m_clipped[m_clipped > 3 ] = 3
  m_clipped[m_clipped < -3] = -3
  #heat.map = heatmap3(m_clipped, Rowv =as.dendrogram(hr), Colv = NA, scale = "none", RowSideColors = mycolhc_c)
  distance_m_c = dist(t(m), method = "euclidean")
  hc <- hclust(distance_m_c, method = "ward.D2")
  if(!is.null(filename)){
    png(filename = filename, pointsize=10, width=1400, height=1400, res=300)
    #svg(filename = filename, width = 15, height = 15)
  }
  if(cluster.col == TRUE){
    col.cluster = as.dendrogram(hc)
  } else{
    col.cluster = NA
  }
  map.col = heatmap3(m_clipped, Rowv =as.dendrogram(hr), Colv = col.cluster, RowSideColors = mycolhc_c, ) 
  if(!is.null(filename)){
    dev.off()
    }
  return(list(
    #heatmap = heat.map,
    heatmap = map.col,
    dendrogram = mycl
  ))
}

#' creates cluster distance mapping. namely the frequency table of every cluster per siRNA
#'
#' @param clustering dendrogram of hierarchical clustering. Output of create clustering method
#'
#' @return cluster frequency table
#' @export
#'
#' @examples
cluster_distance = function(clustering){
  cluster.freq = data.frame(cl = as.vector(clustering), type = names(clustering))
  summ = cluster.freq %>% group_by(cl, type) %>% summarise(count = n())  %>% ungroup() %>% mutate(cl = as.factor(cl))
  return(list(summarised = summ,
              cluster.freq = cluster.freq))
}
#' Creates the cluster representitatives as mean trajectories per cluster.
#'
#' @param cluster.freq cluster frequency table
#' @param features features table used for clustering (as it has the same ordering still)
#' @param dt.data raw time series data
#' @param site.var column identifier of the site.
#' @param track.var column identifier of the track_id
#' @param time.var column identifier of the time variable of the timeseries.
#' @param robust uses median instead of mean, but does not calculate medoids.
#'
#' @return table of average timeseries data.
#' @export
#'
#' @examples
cluster_representatives = function(cluster.freq, features, dt.data, site.var, track.var, time.var, robust = F){
  
  # should still have same ordering as features table
  cluster.freq = mutate(cluster.freq,Image_Metadata_Site = features %>% pull(site.var))
  cluster.freq = mutate(cluster.freq,track_id = features %>% pull(track.var))
  # cluster.freq$Image_Metadata_Site = data.clean$Image_Metadata_Site
  # cluster.freq$objNuc_TrackObjects_Label = data.clean$objNuc_TrackObjects_Label
  
  timeseries.clustered = left_join(cluster.freq, dt.data)
  if(robust){
    summarised_clusters = timeseries.clustered %>%
      mutate(cl = as.factor(cl)) %>%
      group_by_("cl", time.var) %>% summarise(mean_cl_at_t = median(erk.ratio)) 
  }else{
    summarised_clusters = timeseries.clustered %>%
      mutate(cl = as.factor(cl)) %>%
      group_by_("cl", time.var) %>% summarise(mean_cl_at_t = mean(erk.ratio)) 
  }

  return(summarised_clusters)
}
# exp var only nescesary with pooled data, then dt.data and features need a exp identifier to localise median
#' Calculates medoids of each cluster
#'
#' @param dend dendorgram from clustering.
#' @param dt.data raw time series data
#' @param site.var column identifier of the site.
#' @param track.var column identifier of the track_id.
#' @param time.var column identifier of the time variable
#' @param erk.ratio.var column identifier of the erk meausurements
#' @param exp.var column identifier of the site of the experiment to accuratly match trajectories from feature to rawtime series
#' @param features tibble of features per TS
#'
#' @return table containg medoid trajectory per cluster
#' @export
#'
#' @examples
get_medoids = function(dend, features,dt.data, site.var, track.var, time.var, erk.ratio.var, exp.var = NULL){
  
  w = lapply(seq(max(unique(dend))), function(i){
    sums = colSums(as.matrix(dist(features[which(dend == i),])),  na.rm = T)
    a = which(dend == i)[which.min(sums)]
    t = inner_join(dt.data, features[a,]  %>%
                     select(track.var, site.var, exp.var)) %>% 
      ungroup  %>% select(time.var, erk.ratio.var) %>% mutate(cluster = as.factor(i))
    return(t)
  })
  do.call(rbind, w)
}

#'Plots a bar plot of the clutster frequencies. 
#'
#' @param freq.summarised summarised frequencies per siRNA. Output of cluster_distance function
#' @param clustered boolean, indicating if plot should be ordered based on simple clustering ond frequency table.
#'
#' @return ggPlot of cluster frequencies.
#' @export
#'
#' @examples
cluster_dist_plot = function(freq.summarised, clustered = F){
  if(clustered == T){
  freq.summarised %>% spread(cl, count) %>% mutate_if(is.numeric , replace_na, replace = 0) -> temp
  names = temp$type
  m = as.matrix(temp %>% select(-type))
  m = m / rowSums(m)
  rownames(m) = as.character(names)
  distance_m = dist(m, method = "euclidean")
  hr <- hclust(distance_m, method = "ward.D2")
  
  #plot(hr)
  
  as.character(names[hr$order])
  x = seq(length(names[hr$order]))
  names(x) = as.character(names[hr$order])
  ord = dplyr::recode_factor(freq.summarised$type, !!!x)
  freq.summarised$order = as.integer(as.character(ord))
  
  gg = ggplot(freq.summarised,aes(x = reorder(type, order), y = count, fill = cl)) + 
    geom_bar(position = position_fill(reverse = TRUE),stat = "identity") +
    scale_y_continuous(labels = scales::percent_format()) +
    ggplotTheme() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(y = "precentage of trajectories", fill = "Cluster no.", x = "siRNA treatment") 
  }else {
  
  gg = ggplot(freq.summerised,aes(x = type, y = count, fill = cl)) + 
    geom_bar(position = position_fill(reverse = TRUE),stat = "identity") +
    scale_y_continuous(labels = scales::percent_format()) +
    ggplotTheme() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(y = "precentage of trajectories", fill = "Cluster no.", x = "siRNA treatment")
  }
  return(gg)
}

#' Plots the mean trajectories of a clustering
#'
#' @param cluster_mean cluster means as clalculated in the cluster_representative function
#'
#' @return ggPlot of cluster means
#' @export
#'
#' @examples
cluster_mean_plot = function(cluster_mean){
gg = ggplot(cluster_mean, aes(x = Image_Metadata_T, y = mean_cl_at_t, color = cl)) + 
  geom_line() +
  labs(x = "Real Time (min)" , y = "Cytoplasmic to Nuclear intensity ratio ERK", color = "Cluster",  title = "Clustermeans") +
  ggplotTheme()
  return(gg)
}
#' Plots the ouput of the get_medoids function. 
#'
#' @param cluster_median Output of get_medoid function
#'
#' @return ggPlot of cluster medoids
#' @export
#'
#' @examples
cluster_medioid_plot = function(cluster_median){
  gg = ggplot(cluster_median ,aes(x = Image_Metadata_T, y = erk.ratio, color = cluster)) +
    geom_line() +
    labs(x = "Real Time (min)" , y = "Cytoplasmic to Nuclear intensity ratio ERK", color = "Cluster", title = "Medoids") +
    ggplotTheme()
  return(gg)
}

#' composite function creating whole clustering and plots. Saves plots to filename
#'
#' @param features features per trajectories
#' @param data raw timeseries data
#' @param k number of clusters
#' @param filename filname where to save images 
#' @param stim.var column identifier of Stimulation_treatment
#' @param site.var column identifier of site
#' @param track.var column identifier of track_id
#' @param time.var column identifier of time variable
#' @param remove.var list of features to remove
#' @param exp.var column identifier of the experiment identifier
#' @param erk.ratio.var column identifier of erk measurement
#' @param normalize_plot boolean, normatlizes baseline if T
#'
#' @return None.
#' @export
#'
#' @examples
cluster = function(features, data, k, filename, stim.var, site.var, track.var, time.var, remove.var, exp.var, erk.ratio.var, normalize_plot = F){

  clustering = create_clustering(features, stim.var,site.var, track.var, time.var, k, remove.var = remove.var, filename =paste0(filename, ".png"))
  freq = cluster_distance(clustering$dendrogram)
  if(normalize_plot){
    dt.data.m = data %>% filter(get(time.var) < 5) %>% group_by(.dots = c(site.var, track.var)) %>% summarise(diff = 1 - mean(get(erk.ratio.var)))
    data = left_join(data, dt.data.m) %>% mutate(erk.ratio.n = get(erk.ratio.var) + diff) %>% select(-erk.ratio.var) %>%
      rename_(.dots = setNames("erk.ratio.n", erk.ratio.var))
  }
  cluster_mean = cluster_representatives(cluster.freq = freq$cluster.freq, features = features, dt.data = data, site.var, track.var, time.var)
  cluster_median = get_medoids(clustering$dendrogram, features, data, site.var, track.var,time.var, erk.ratio.var, exp.var)
  
  gg = cluster_dist_plot(freq$summarised, clustered = T)
  ggsave(paste0(filename,"_dist" ,".svg"), gg + labs_pubr(base_size = 20), device = "svg", width = 12, height = 12)
  gg = cluster_mean_plot(cluster_mean)
  ggsave(paste0(filename, "_mean",".svg"), gg, device = "svg", width = 12, height = 12)
    gg = cluster_medioid_plot(cluster_median)
  ggsave(paste0(filename, "_medoid",".svg"), gg, device = "svg", width = 12, height = 12)
}

# Script not function, gets not called if sourced.
if(FALSE){
  
  ##!!EXAMPLE
  load("data/features.three.n.RData")
  dt.data = load("data/data.three.RData")
  source("data.R")
  features.clean = features.three.n
  dt.data = data.three
  
  features.clean.pooled = features.clean %>% pool_treatments(., stim.var)

  
  names = features.pooled.n %>% pull(var = stim.var)
  cluster.data.n = features.pooled.n %>% ungroup()  %>% 
    select(-site.var, -track.var) %>% #, -dip.time.min, -mean, -sd, -FMWH.left, -FMWH.right, -noise) %>% 
    select_if(is.numeric) 
  # m = as.matrix(cluster.data.n)
  # rownames(m) = names
  k = 5
  
  clustering = create_clustering(data = features.pooled.n, stim.var, site.var, track.var, time.var, k)
  freq = cluster_distance(clustering$dendrogram)
  cluster_mean = cluster_representatives(cluster.freq = freq$cluster.freq, features = features.clean, dt.data = dt.data, site.var, track.var, time.var)
  cluster_median = get_medoids(clustering$dendrogram, features.pooled.n, dt.data, site.var, track.var,time.var)
  
  cluster(features.pooled.n, dt.data, k, "filename", stim.var, site.var, track.var, time.var, remove.var, exp.var, erk.ratio.var, normalize_plot = F)
}