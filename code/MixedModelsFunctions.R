
if (!all(c("sommer", "lme4")%in% installed.packages())) {
  if(!"sommer" %in% installed.packages()) {install.packages("sommer")
    library(sommer)}
  if(!"lme4" %in% installed.packages()) {install.packages("lme4")
    library(lme4)}
} else {library(sommer); library(lme4)}


##### #### Analise de modelos mistos para ensaio com delineamento de blocos aumentados
#### e completos - sommer package
analyzeTrial.sommer <- function(x){
  if(any(x$studyDesign == "DBA")){
    #### Augmented Blocks
    modfit <- mmer(y ~ check + rep,
                   random = ~ clone:new,
                   data = x,
                   tolparinv = 1e-6,
                   verbose = F,
                   getPEV = T)
  } else {
    #### Complete Block Design
    modfit <- mmer(y ~ rep,
                   random = ~ clone,
                   data=x,
                   tolparinv = 1e-6,
                   verbose = F,
                   getPEV = T)
  }
  return(modfit)
}


getVarComp.sommer <- function(model){
  unlist(model$sigma) %>%
    data.frame(effect = c("Clone", "Residual"), VarEstimate = .) %>%
    spread(key = effect, value = VarEstimate, fill=NA)
}
