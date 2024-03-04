library(hht)#CEEMD
library(ggplot2)


hanr_remd <- function() {
  obj <- harbinger()
  obj$sw_size <- NULL
  
  class(obj) <- append("hanr_remd", class(obj))
  return(obj)
}

fit.hanr_remd <- function(obj, serie, noise = 0.1, ...) {
  if(is.null(serie)) stop("No data was provided for computation",call. = FALSE)
  
  serie <- stats::na.omit(serie)
  w=30
  # original era 0.001
  noise.amp = noise
  trials=5
  serie <- ts(serie)
  
  id <-  1:length(serie)
  obj$sw_size <-  length(serie)
  ceemd.result <- CEEMD(serie, id, noise, trials)
  
  obj$model <- ceemd.result
  
  return(obj)
}


detect.hanr_remd <- function(obj, serie, ...) {
  if(is.null(serie)) stop("No data was provided for computation",call. = FALSE)
  
    vec <- vector()
  for (n in 1:obj$model$nimf){
    vec[n] <- fc_rugosidade(obj$model[["imf"]][,n])
  }
  ##Curvatura máxima
  res <- transform(fit_curvature_max(), vec)
  div <- res$x

  soma_high_freq <- obj$model[["imf"]][,1]
  
  for (k in 2:div){
    soma_high_freq <- soma_high_freq + obj$model[["imf"]][,k]}
  diff_soma <- c(NA, diff(soma_high_freq)) 
  
  media_test <- median(abs(diff_soma), na.rm= TRUE)
  distancia <- abs(diff_soma) -  media_test
  outliers_emc <- which(abs(distancia)>2.698*sd(distancia, na.rm=TRUE))

  ts <- ts_data(soma_high_freq, 0)
  io <- ts_projection(ts)
  model <- ts_arima()
  model <- fit(model, x=io$input, y=io$output)
  adjust <- predict(model, io$input)
  adjust <- as.vector(adjust)
  
  diferencas_absolutas <- abs(adjust - soma_high_freq)
  probabilidades <-(1 - diferencas_absolutas / max(diferencas_absolutas))
  outliers_arima <- which(abs(probabilidades)<2.698*sd(probabilidades, na.rm=TRUE))

  intersecao <- intersect(outliers_emc, outliers_arima)
  
  obj$anomalies[1:obj$sw_size] <- FALSE
  
  if (!is.null(intersecao) & length(intersecao) > 0) {
    obj$anomalies[intersecao] <- TRUE
  }
  
  detection <- obj$har_restore_refs(obj, anomalies = obj$anomalies)
  
  return(detection)
  
  
}

##Função de rugosidade
fc_rugosidade <- function(x){
  firstD = diff(x)
  normFirstD = (firstD - mean(firstD)) / sd(firstD)
  roughness = (diff(normFirstD) ** 2) / 4
  return(mean(roughness))
}