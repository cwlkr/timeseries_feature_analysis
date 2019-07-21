#' @author Cédric Walker
#' !!!Do Not RUN
#' Script used to create annotation for random forest and train random forest on the annotations
#' See thesis Chapter Identification of introduced error in feature calculations 
#' 

# create oulier detection decosion tree
require(TSexplore)
require(tidyverse)
require(rpart.plot)
require(rpart)

source("plot_features_per_series.R")
#cyto.erk = "objCyto_Intensity_MeanIntensity_imErk"
time.var = "RealTime"
stim.var = "Stimulation_treatment"
stim.time.var = "Stimulation_time"
group.var = c("Image_Metadata_Site", "objNuc_TrackObjects_Label")
erk.ratio.var = "erk.ratio"
# manually annotated set.
## for annotation-------------
#training was done on the default 
source("data.R")
# dt.data = load_data(experiment  = "/Coralie/NIH3T3/siPOOLs/20190311_systII_siPOOLs_plate1_and_2_singlePulse/20190311_102308_731")$data
# fwrite(dt.data, "bs_forest/20190311_systII_siPOOLs_plate1_and_2_singlePulse_random_forest_annotation_prepare.csv")
fread("bs_forest/20190311_systII_siPOOLs_plate1_and_2_singlePulse_random_forest_annotation_prepare.csv")
data = dt.data %>% unite("unique" , c(track.var, Image_Metadata_Site), sep = "_") %>% filter(!(Stimulation_treatment %like% "+CTRL"))

set.seed(156)
k = sample(unique(data$unique), 400)
fwrite(as.list(k) , "kconfic20190311_plate1and2.csv")
pdf("bs_annotator_information.pdf20190311_plate1and2.csv")

for(i in 1:length(k)){
   single_timeseries_coralie = data %>%  filter(unique == k[i]) %>%mutate(nuc_erk = objCyto_ring_Intensity_MeanIntensity_imErk/objNuc_Intensity_MeanIntensity_imErk) %>% select(nuc_erk) %>% melt() %>% pull()
   
   baseline = series_baseline(single_timeseries_coralie, 10, robust = T)
   dip_amp = dip_amplitude(single_timeseries_coralie, stim.time = 10, baseline)
   max_amp = max_amp_wrapper(single_timeseries_coralie, stim.time =  stim.times, baseline = baseline, dip.amp = dip_amp)
   fwhm = fwhm_wrapper(single_timeseries_coralie,  basal = baseline$before, dip_amp)
   dec = decay_half_max(single_timeseries_coralie, max.amp = max_amp, fwhm = fwhm)
   growth = growth_half_max(single_timeseries_coralie, max.amp = max_amp, fwhm = fwhm)
   slope.after.max = afterpeakdecay(single_timeseries_coralie, max_amp, 3)
                                     
   series_plot_feat(single_timeseries_coralie, baseline, dip_amp, max_amp, fwhm, dec, growth, slope.after.max, stim.times)
   single_timeseries_coralie_treat = data %>% filter(unique == k[i]) %>%
     select(Stimulation_treatment) %>% pull()
   title(paste(as.character(single_timeseries_coralie_treat)[1], k[i] , sep = "_"))
   
 }
 dev.off()

 annot = fread(file = "bs_forest/annot_20190311_plate1and2_classes.csv")
 data.annot = data %>% filter(unique %in% k)
 data.annot$id = data.annot$unique
 data.annot = data.annot %>% separate(unique, into = c('track_id', 'Image_Metadata_Site'), sep = "_") %>% mutate(track_id = as.numeric(track_id),
                                                                                                                Image_Metadata_Site = as.numeric(Image_Metadata_Site)) 
 # make shure they are in the  right order
 features = all_features_df(data.annot %>% filter(id  == k[1])  %>% select(erk.ratio) %>% melt() %>% pull(),stim.time = 10)
 for(i in k[-1]){
   series = data.annot %>% filter(id  == i)  %>% select(erk.ratio) %>% melt() %>% pull()
   features  =  rbind(features, all_features_df(series, 10))
 }
 
 f_data = features #%>% mutate_if(is.numeric,function(.){(. - mean(., na.rm = T))/sd(., na.rm = T)})
 f_data[is.na(f_data)] = 1000
 f_data$class = as.factor(annot$annot_class +1) 
 f_data$Stimulation_treatment = NULL

testind = sample(1:nrow(features), 300, replace = F)
test = f_data
test = test[-testind,]
testclasses = test$class
test$class = NULL
train = f_data[testind,]
train$Image_Metadata_Site = NULL
train$objNuc_TrackObjects_Label = NULL


 require(randomForest)
 fit <- randomForest(class ~ ., data = train,ntree = 10000 )
 pred.rnd = predict(fit, test)
 table(pred = pred.rnd, true = testclasses)
 mean(pred.rnd == testclasses)
 print(fit) # view results
 importance(fit) # importance of each predictor 
 
 
 #model works, train it with both

 fit <- randomForest(class ~ ., data = f_data,ntree = 10000 )
 pred.rnd = predict(fit, test)
 table(pred = pred.rnd, true = testclasses)
 mean(pred.rnd == testclasses)
 print(fit) # view results
 importance(fit) # importance of each predictor 

 # safe model !
 write_rds(fit, "bs_forest/bs_forest.classifier_data20190311_plate1and2")
 