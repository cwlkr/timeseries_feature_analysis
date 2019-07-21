#'  @Author Cédric Walker



#' Title
#'
#' @param data 
#' @param ci.lvl 
#' @param stimulus.rug 
#' @param stim.var 
#' @param erk.ratio.var 
#' @param time.var 
#' @param stim.color 
#' @param alpha 
#' @param xlab 
#' @param ylabel 
#' @param vlines 
#'
#' @return
#' @export
#'
#' @examples
create_plot = function(data, ci.lvl, stimulus.rug, stim.var, erk.ratio.var, time.var, stim.color = "red",  alpha = 0.1, xlab = "Real Time (min)", ylabel = "ERK-KTR Cytoplasmic ot Nuclear Ratio", vlines = F){
  source("data.R")
  a = ci.lvl
  #create unique id for each track..
  #data[,unique := paste(Image_Metadata_Site, objNuc_TrackObjects_Label, sep = "_")]
  #calc erk ratio
  # for each grouping calc mean, min, max, ci upper lower
  data.summary <- data %>% ungroup() %>% group_by(.dots = c(stim.var, time.var)) %>% summarise(ymin_ctn = min(get(erk.ratio.var)) , ymax_ctn = max(get(erk.ratio.var)), mean_ctn = mean(get(erk.ratio.var)),
          lower = mean(get(erk.ratio.var)) - (qnorm(1-(a/2)) * sd(get(erk.ratio.var))/sqrt(n())) , upper =  mean(get(erk.ratio.var)) + (qnorm(1-(a/2)) * sd(get(erk.ratio.var))/sqrt(n())))
    #get also rug of pulses
  
  #check if CTRL and treat differently and maybe check if everything is ctrl.
  ctrl = data.summary %>% filter(get(stim.var) %like% "CTRL")
  #ctrl.neg = data.summary %>% filter(get(stim.var) %like% "\\-CTRL")
  #ctrl.pos = data.summary %>% filter(get(stim.var) %like% "\\+CTRL")
  #ctrl = ctrl %>% ungroup %>% mutate(Stimulation_treatment =  forcats::fct_rev(as.factor(get(stim.var))))
  
  if(nrow(ctrl) == nrow(data.summary)){
    data.plot = ctrl
  }else{
    data.plot = data.summary %>% filter(!(get(stim.var) %like% "CTRL"))
  }
  
  ggp = ggplot(data.plot, aes(x = get(time.var), y = mean_ctn, group = get(stim.var)))+ 
    geom_ribbon(alpha = alpha, mapping = aes(ymin = lower, ymax = upper)) +
    geom_line(aes(color = get(stim.var)), size = 1.2) + 
    labs(x = time.var, y = ylabel, legend = "stim.var")  + ggplotTheme() + theme(legend.title = element_blank()) +
    geom_rug(data = data.frame(stim.times = stim.times), aes(x=stim.times,y =NULL, group = NULL), color = stim.color)
    
    if(!(nrow(ctrl) == nrow(data.summary))){
      ggp = ggp + geom_line(data = ctrl, aes( group = get(stim.var), linetype=get(stim.var)),  color = "black",  size = 1.2) +
      #geom_line(data = ctrl.neg, aes( group = get(stim.var)), linetype="dashed", color = "black",  size = 1.2) +
      #scale_linetype_manual(values = c("pos", "neg")) +
      geom_ribbon(data = ctrl,alpha = alpha, mapping = aes(ymin = lower, ymax = upper)) + 
      scale_color_brewer(palette = "Set1")
    }
    
    if(vlines){
      ggp = ggp + geom_vline(data = data.frame(stim.times = stim.times), aes(xintercept=stim.times, group = NULL), color = stim.color, linetype="dotted")
    }
    
    return(ggp)
  }

custom_rug = function(ggplot.obj, stim.times, stim.color){
  rug = ggplot.obj + geom_rug(aes(x=stim.times[1],y = NULL), color = stim.color)
  for (ta in stim.times[-1]){
    rug = rug + geom_rug(aes(x=ta,y = NULL), color = stim.color)
  }
  return(rug)
}


#create_plot(data, 0.05, 10, nuc.erk = nuc.erk, cyto.erk = cyto.erk, time.var = time.var, stim.var = stim.var)

#' Title
#'
#' @param dt.data 
#' @param pdf.filename 
#' @param meta.grouping 
#' @param stim.times 
#' @param time.var 
#' @param stim.var 
#' @param erk.ratio.var 
#' @param vlines 
#'
#' @return
#' @export
#'
#' @examples
plot_average_per_group <- function(dt.data, pdf.filename, meta.grouping, stim.times, time.var, stim.var, erk.ratio.var, vlines) {
  pdf(file = pdf.filename)
  ngroups = dt.data %>% select(meta.grouping)  %>% summarise(n = length(unique(get(meta.grouping)))) %>% pull()
  for(i in (0:(ngroups-1))){
    
    #stim_vec = as.numeric(pull(stim_vec))
    gg = create_plot(data = setDT(dt.data)[get(meta.grouping) == i],
                     ci.lvl = 0.05,
                     stimulus.rug = stim.times,
                     time.var = time.var,
                     stim.var = stim.var,
                     erk.ratio.var = erk.ratio.var,
                     vlines = vlines)
    plot(gg)
  }
  gg = create_plot(data = setDT(dt.data)[get(stim.var) %like% "CTRL"],
                   ci.lvl = 0.05,
                   stimulus.rug = stim.times,
                   time.var = time.var,
                   stim.var = stim.var,
                   erk.ratio.var = erk.ratio.var,
                   vlines = vlines)
  plot(gg)
  dev.off()
}


