#'  @Author Cédric Walker


#' Title
#'
#' @param data 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
normalize_frame = function(data, ...){
  require(tidyverse)
  n = ncol(data)
  #data %>% ungroup  %>% mutate_if(is.numeric, .funs = function(.) (mean(., na.rm=T) - .) /sd(., na.rm=T))
  data %>% ungroup  %>% mutate_at(.vars = vars(...), .funs = function(.) if(is.numeric(.)){(. - mean(., na.rm=T)) /sd(., na.rm=T)}else{.})
}


