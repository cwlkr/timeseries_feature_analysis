#' @author  Cédric Walker
#' see thesis chapter feature importance and selection 


# calculate feature importance with a random forest
library(randomForest)
library(ggpubr)
#rnd.features = data.clean
load("data/features.three.n.RData")
source("load_multiple.R")
features.clean = pool_treatments(features = features.three.n, stim.var)
rnd.features = features.clean
# create own class for NA data by setting it to 1000
rnd.features[is.na(rnd.features)] = 1000
rnd.features = select(rnd.features, -site.var)
rnd.features = select(rnd.features, -track.var, -amplitude.time.max, -sd, -noise, -exp, -mean)
rnd.features = mutate(rnd.features,Stimulation_treatment = as.factor(get(stim.var)))
forest_class = randomForest(Stimulation_treatment~.,rnd.features, importance =T)
coeff = as.data.table(importance(forest_class))
coeff$coeff = rownames(importance(forest_class))


c = ggplot(coeff, aes(y = MeanDecreaseGini, x = coeff, fill = MeanDecreaseGini)) + 
  geom_bar(stat = "identity") + theme_pubr(legend = "right")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Features")

ggsave(c, filename = "plotsforthesis/final/forest.svg", device = "svg")


# Load and run Boruta
# run boruta
names = features.clean$Stimulation_treatment
boruta.data = features.clean %>% ungroup()  %>% 
  select(-Image_Metadata_Site, -track_id, -noise, -baseline.after.slope, -sd, -mean, -exp) %>% mutate(Stimulation_treatment = as.factor(Stimulation_treatment)) %>% 
  replace(., is.na(.), 10000)

b = Boruta::Boruta(Stimulation_treatment~., data = boruta.data)

