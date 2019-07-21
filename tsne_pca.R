require(Rtsne)
source("data.R")
load("data/features.three.n.RData")

data.pooled.n = features.three.n %>% mutate(Stimulation_treatment = case_when(  Stimulation_treatment %like% "\\+CTRL" ~ "+CTRL",
  Stimulation_treatment %like% "\\-CTRL" ~ "-CTRL",
  Stimulation_treatment %like% "WT" ~ "WT",
  TRUE ~ as.character(Stimulation_treatment)))

# try tsne stuff
labels = na.omit(data.pooled.n) %>% pull(stim.var)
train<- data.pooled.n %>% select(-site.var, -track.var, -stim.var)## Choose the train.csv file downloaded from the link above
## Curating the database for analysis with both t-SNE and PCA

labels<-as.factor(labels)
## for plotting
colors = rainbow(length(unique(labels)))
names(colors) = unique(labels)

## Executing the tsne algorithm on curated data
tsne <- Rtsne(na.omit(train), dims = 2, perplexity=50, verbose=TRUE, max_iter = 2500, eta = 300, num_threads = 4)

## Plotting tsne
x = data.frame(tsne$Y)
x$type = labels
ggplot(x, aes(x = X1, y = X2, color = type)) +
  geom_point() + theme_classic()


#pca
d = data.pooled.n[sample(nrow(data.pooled.n),5000),]
names = d$Stimulation_treatment
cluster.data = d %>% ungroup()  %>% 
  select(-Image_Metadata_Site, -track_id, -noise, -baseline.after.slope, -sd, -mean, -exp) %>% select_if(is.numeric) 
cluster.data.n = normalize_frame(cluster.data)
m = as.matrix(cluster.data.n)

m.pca = prcomp(na.omit(m),center = TRUE)
x = ggbiplot::ggbiplot(m.pca, groups=na.omit(d)$Stimulation_treatment, ellipse = T, alpha = 1)
x + scale_y_continuous(limits = c(-5,5)) + scale_x_continuous(limits = c(-5,5))

ggsave(x, filename = "plotsforthesis/pca.svg", device = "svg", width = 12, height = 12)

