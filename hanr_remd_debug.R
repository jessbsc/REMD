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
  
  ######## plot todas as decomposicoes. 
  # PODE APAGAR SE NAO QUISER FAZER O PLOT
  # num_colunas_ifms <- obj$model$nimf
  # colunas_para_adicionar <- list()
  # for (i in 1:num_colunas_ifms) {
  #   nome_coluna <- paste(i, "° imf")
  #   colunas_para_adicionar[[nome_coluna]] <- obj$model[["imf"]][,i]
  # }
  # colunas_para_adicionar[['residue']] <- obj$model[["residue"]]
  # # Usar ts.union() com todas as colunas de 'ifms'
  # resultado_union <- do.call(ts.union, colunas_para_adicionar)
  # # Plotar o resultado
  # if (num_colunas_ifms < 10){
  #   plot(resultado_union, main = "Decomposição Completa")
  # }else{
  #   plot(resultado_union[,0:10], main = "Decomposição Completa")
  # }

  
  ##calculando rugosidade pra cada imf
  vec <- vector()
  for (n in 1:obj$model$nimf){
    vec[n] <- fc_rugosidade(obj$model[["imf"]][,n])
  }
  ##Curvatura máxima
  res <- transform(fit_curvature_max(), vec)
  div <- res$x
  #plot(vec, t="p", col=ifelse(vec==vec[div], 'red', 'blue'))
  #lines(vec, type="b", pch=NA, lty=3, col="black")
  soma_high_freq <- obj$model[["imf"]][,1]
  
  for (k in 2:div){
    soma_high_freq <- soma_high_freq + obj$model[["imf"]][,k]}
  #plot(soma_high_freq, type="l")
  #title(main = "Soma das maiores frequências")
  # 
  # identificacao dos pontos de anomalia
  diff_soma <- c(NA, diff(soma_high_freq)) 
  # plot(diff_soma, type="l")
  # title(main = "Diferencial da soma das mariores frequencias - variacao entre os elementos consecutivos")
  
  media_test <- median(abs(diff_soma), na.rm= TRUE)
  distancia <- abs(diff_soma) -  media_test
  outliers_emc <- which(abs(distancia)>2.698*sd(distancia, na.rm=TRUE))
  
  
  # Arima Filter 
  ts <- ts_data(soma_high_freq, 0)
  io <- ts_projection(ts)
  model <- ts_arima()
  model <- fit(model, x=io$input, y=io$output)
  adjust <- predict(model, io$input)
  adjust <- as.vector(adjust)
  
  # Calculo da probabilidade inversa ( ou seja, quem esta mais proximo de 1 tem mais chance de ser outlier - desse jeito é melhor para calcular o outliers abaixo)
  diferencas_absolutas <- abs(adjust - soma_high_freq)
  probabilidades <-(1 - diferencas_absolutas / max(diferencas_absolutas))
  #print(probabilidades) 
  outliers_arima <- which(abs(probabilidades)<2.698*sd(probabilidades, na.rm=TRUE))
  
  #print(outliers_arima)
  
  # Arima plot 
  # grf <- har_plot(model, soma_high_freq, detection, as.logical(soma_high_freq))
  # # <- grf + geom_vline(xintercept = 75, col = "black", linetype = "dashed")
  # grf <- grf + geom_line(aes(y=adjust), linetype = "dashed", col="darkblue")
  # grf <- grf + geom_point(aes(y=adjust), size=0.25, col="darkblue")
  # grf <- grf
  # plot(grf)
  # 
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