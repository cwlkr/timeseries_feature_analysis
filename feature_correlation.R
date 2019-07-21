# plot feature correlation matrix
# corr ---------------
source("data.R")
source("normalize_data_frame.R")
source("load_multiple.R")
load("data/features.three.n.RData")

exp.var = "exp"
site.var = "Image_Metadata_Site"
track.var = "track_id"
time.var = "Image_Metadata_T"
stim.var = "Stimulation_treatment"
erk.ratio.var = "erk.ratio" 
features.clean.pooled = pool_treatments(features.three.n, stim.var) %>% filter(!get(stim.var) %like% "\\+CTRL")


plot_feature_correlation(features.clean.pooled, track.var, site.var, exp.var, T)
plot_feature_correlation <- function(features.clean.pooled, track.var, site.var, exp.var, show_labels = F, filename = NULL) {
  
  features.clean.pooled %>% ungroup %>% select_if(is.numeric) %>% select(-track.var, -site.var, -exp.var) %>% na.omit %>%
    select(-sd, -mean, -baseline.after.slope, -noise, - amplitude.time.max)  %>% cor %>% data.frame()-> cor.data 
  cor.data$type = rownames(cor.data)
  cor.long = gather(cor.data, "type1","Correlation",-type) %>% mutate(Correlation = round(Correlation,2))
  cor_p = ggplot(cor.long, aes(y = reorder(type, Correlation), x = reorder(type1, Correlation), fill = Correlation)) + 
    geom_tile() + scale_fill_gradient2(low = "blue", high= "red", mid = "white") + theme(text =element_text(size = 20)) +labs_pubr(23) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Correlation of Features" , x = "" , y = "")
  if(show_labels){
    cor_p = cor_p +  geom_text(aes(label = Correlation), color = "black", size =   4) 
  }
  
  plot(cor_p)
  if(!is.null(filename)){
    ggsave(cor_p, filename = "plotsforthesis/final/correlation_matrix.svg", device = "svg", width = 14, height = 12, units = "in")
    
  }
}
