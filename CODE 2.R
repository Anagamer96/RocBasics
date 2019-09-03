###############################################################################
# Funcion auxiliar para hacer un gráfico de las distintas poblaciones

poblaciones1 <- function(casos, controles){
  
  #############################################################################
  # casos: Valores de los individuos enfermos
  # Controles: Valores de los individuos sanos
  #############################################################################

  if(!is.numeric(casos))     stop("Los casos deben ser numericos")
  if(!is.numeric(controles)) stop("Los controles deben ser numericos") 
  
  par(mfrow=c(2,1))
  hist(casos, col = "red",
       xlim = c(min(c(controles, casos)),max(c(controles, casos))),
       main = "Casos", xlab = "Valor del biomarcador", ylab = "Frecuencia")
  hist(controles, col = "blue",
       xlim = c(min(c(controles, casos)),max(c(controles, casos))),
       main = "Controles", xlab = "Valor del biomarcador", ylab = "Frecuencia")
  par(mfrow=c(1,1))
}

set.seed(285899)
poblaciones1( casos = rnorm(50, 7, 2), controles = rnorm(50, 4, 2))



###############################################################################
# Funcion auxiliar para hacer un gráfico de las distintas poblaciones

poblaciones2 <- function(casos, controles){
  
  #############################################################################
  # casos: Valores de los individuos enfermos
  # Controles: Valores de los individuos sanos
  #############################################################################
  
  if(!is.numeric(casos))     stop("Los casos deben ser numericos")
  if(!is.numeric(controles)) stop("Los controles deben ser numericos") 
  
  plot(density(casos), col = "red",
       xlim = c(min(c(controles, casos)),max(c(controles, casos))), 
       ylim = c(0,1.05*max(c(max(density(casos)$y)),max(density(controles)$y))),
       main = "Poblaciones", xlab = "Valor del biomarcador", ylab = "")
  lines(density(controles), col = "blue")
}

set.seed(285899)
poblaciones2(casos = rnorm(50, 7, 2), controles = rnorm(50, 4, 2))



###############################################################################
# Funcion auxiliar para el calculo de la sensibilidad y especificidad 

calculo_SeSp <- function (casos, controles){
  
  #############################################################################
  # casos: Valores de los individuos enfermos
  # Controles: Valores de los individuos sanos
  #############################################################################
  
  if(!is.numeric(casos))     stop("Los casos deben ser numericos")
  if(!is.numeric(controles)) stop("Los controles deben ser numericos") 
  
  # Creacion de variables auxiliares y manipulacion de argumentos
  enf <- length(casos)
  san <- length(controles)
  casos <- sort(casos)
  controles <- sort(controles)
  cc <- sort(unique(c(casos, controles)))
  aux1 <- min(cc) - (0.01*min(cc))^2
  aux2 <- max(cc) + 0.01*max(cc)
  cc <- c(aux1, cc, aux2)
  
  roc <- data.frame(cc)
  # Cálculo de sensibilidad y especificidad para diferentes puntos de corte
  for(i in 1:length(cc)){
    roc$Se[i] <- sum(casos >= roc$cc[i])/enf
    roc$Sp[i] <- sum(controles < roc$cc[i])/san
    # roc$Se[i] <- sum(casos > roc$cc[i])/enf
    # roc$Sp[i] <- sum(controles < roc$cc[i])/san
  }
  
  return(roc)
}

set.seed(285899)
calculo_SeSp(casos = rnorm(50, 7, 2), controles = rnorm(50, 4, 2))



###############################################################################
# Funcion auxiliar para el calculo de la sensibilidad y especificidad tomando 
# la distribucion suavizada de los datos recogidos

calculo_SeSp_suavizado <- function (casos, controles){
  
  #############################################################################
  # casos: Valores de los individuos enfermos
  # Controles: Valores de los individuos sanos
  #############################################################################
  
  if(!is.numeric(casos))     stop("Los casos deben ser numericos")
  if(!is.numeric(controles)) stop("Los controles deben ser numericos") 
  
  # Creacion de variables auxiliares y manipulacion de argumentos
  dens_casos <- density(casos)
  dens_control <- density(controles)
  dens_casos$y <- dens_casos$y/sum(dens_casos$y)
  dens_control$y <- dens_control$y/sum(dens_control$y)
  cc <- sort(unique(c(dens_casos$x, dens_control$x)))
  aux1 <- min(cc) - (0.01*min(cc))^2
  aux2 <- max(cc) + 0.01*max(cc)
  cc <- c(aux1, cc, aux2)
  roc <- data.frame(cc)
  
  # Cálculo de sensibilidad y especificidad para diferentes puntos de corte
  for(i in 1:length(cc)){
    roc$Se[i] <- sum(dens_casos$y[dens_casos$x >= roc$cc[i]] )
    roc$Sp[i] <- sum(dens_control$y[dens_control$x < roc$cc[i]])
  }
  
  return(roc)
}

set.seed(285899)
calculo_SeSp_suavizado(casos = rnorm(50, 7, 2), controles = rnorm(50, 4, 2))



###############################################################################
# Funcion auxiliar para el calculo del AUC

calculo_auc <- function(roc){
  
  #############################################################################
  # roc: objeto creado por la funcion calculo_SeSp, calculo_SeSp_suavizado o
  # calculo_SeSp_parametrico. Contiene datos de la sensibilidad y especificidad
  #############################################################################
  
  if(!is.data.frame(roc)) stop("Parametro roc no valido")

  roc <- roc[nrow(roc):1,]
  auc <- 0
  for(i in 2:(nrow(roc))){
    auc <- auc + ((1-roc$Sp[i])-(1-roc$Sp[i-1]))*roc$Se[i-1]+((1-roc$Sp[i])-(1-roc$Sp[i-1]))*(roc$Se[i]-roc$Se[i-1])/2
  }
  
  return(auc)
}

set.seed(285899)
calculo_auc(calculo_SeSp(casos = rnorm(50, 7, 2), controles = rnorm(50, 4, 2))) 
calculo_auc(calculo_SeSp_suavizado(casos = rnorm(50, 7, 2), controles = rnorm(50, 4, 2))) 



###############################################################################
# Funcion auxiliar para el calculo del intervalo de confianza del AUC 
# utilizando bootstrap

ICauc_bootstrap <- function(casos, controles, Nsim, alpha){
  
  #############################################################################
  # casos: datos recogidos de los enfermos
  # controles: datos recogidos de los sanos
  # Nsim: numero de simulaciones que se van a utilizar en bootstrap
  # alpha: nivel de confianza
  #############################################################################
  
  if(!is.numeric(casos))     stop("Los casos deben ser numericos")
  if(!is.numeric(controles)) stop("Los controles deben ser numericos") 
  if(!is.numeric(Nsim))      stop("Los casos deben ser numericos")
  if(alpha<0 | alpha>1)      stop("Los controles deben ser numericos") 
  if(Nsim > 10000)           warning("Si Nsim es muy grande puede tardar mucho")
  
  enf <- length(casos)
  san <- length(controles)
  controles_sim <- matrix(0, ncol = Nsim, nrow = san)
  casos_sim <- matrix(0, ncol = Nsim, nrow = enf)
  for (i in 1:Nsim) {
    controles_sim[,i] <- sample(controles, san, replace = TRUE)
    casos_sim[,i] <-  sample(casos, enf, replace = TRUE)
  }
  auci <- rep(0,Nsim)
  for (i in 1:Nsim){
    roci <- calculo_SeSp(casos_sim[,i], controles_sim[,i])
    auci[i] <- calculo_auc(roci)
  }
  return(quantile(auci, c(alpha/2, 1-alpha/2))) # PFFF
}

set.seed(285899)
ICauc_bootstrap(casos = rnorm(50, 7, 2), controles = rnorm(50, 4, 2), Nsim = 30, alpha = 0.05)



###############################################################################
# Funcion auxiliar para el calculo del pAUC

calculo_pauc <- function(roc, AUC){
  
  #############################################################################
  # roc: objeto creado por la funcion calculo_SeSp contiene datos de la 
  #      sensibilidad y especificidad
  # AUC: rango de valores del complementario de la especificidad en el que se 
  #      quiere calcular el area
  #############################################################################
  
  if(!is.numeric(AUC) | length(AUC)!=2) 
    {stop("AUC debe ser un vector de longitud 2")}
  if(AUC[1]<0 | AUC[1]>1 | AUC[2]<0 | AUC[2]>1) 
    {stop("AUC no valido. Debe ser mayor que 0 y menor que 1")}
  if(!is.data.frame(roc) | ncol(roc)<3) 
    {stop("Argumento roc no valido")} 
  
  roc <- roc[nrow(roc):1,]
  pauc <- 0
  indices <- which( (1-roc$Sp)>AUC[1] & (1-roc$Sp)<AUC[2] )
  for(i in (min(indices)+1):(max(indices))){
    pauc <- pauc + ((1-roc$Sp[i])-(1-roc$Sp[i-1]))*roc$Se[i-1]+((1-roc$Sp[i])-(1-roc$Sp[i-1]))*(roc$Se[i]-roc$Se[i-1])/2
  }
  pauc <- pauc /(AUC[2]-AUC[1])  

    return(pauc)
}

set.seed(285899)
calculo_pauc( calculo_SeSp(casos = rnorm(50, 7, 2), controles = rnorm(50, 4, 2)), AUC = c(0.5,0.8)) 
calculo_pauc( calculo_SeSp_suavizado(casos = rnorm(50, 7, 2), controles = rnorm(50, 4, 2)), AUC = c(0.5,0.8)) 


###############################################################################
# REPASAR
# Funcion auxiliar para calcular los intervalos de confianza de la sensiblidad
# y la especificidad

IC_SeSp <- function(roc, enf, san, alpha){
  
  #############################################################################
  # roc: objeto creado por la funcion calculo_SeSp contiene datos de la 
  #      sensibilidad y especificidad
  #############################################################################
  
  if(!is.data.frame(roc) | ncol(roc)<3) stop("Argumento roc no valido") 
  if(!is.numeric(enf))                  stop("enf debe ser un numerico")
  if(!is.numeric(san))                  stop("san debe ser un numerico")
  if(alpha<0 | alpha>1) stop("alpha no valido. Debe ser mayor que 0 y menor que 1")
  
  for(i in 1:nrow(roc)){
    if (roc$Se[i] != 0 & roc$Se[i] != 1){
      roc$ciiSe[i] <- ICnorm(roc$Se[i]*enf, enf, alpha)[2]
      roc$cisSe[i] <- ICnorm(roc$Se[i]*enf, enf, alpha)[3]
    }else{
      roc$ciiSe[i] <- roc$Se[i]  
      roc$cisSe[i] <- roc$Se[i] }
    if (roc$Sp[i] != 0 & roc$Sp[i] != 1){
      roc$ciiSp[i] <- ICnorm(roc$Sp[i]*san, san, alpha)[2]
      roc$cisSp[i] <- ICnorm(roc$Sp[i]*san, san, alpha)[3]
    }else{
      roc$ciiSp[i] <- roc$Sp[i]  
      roc$cisSp[i] <- roc$Sp[i] }
    if (roc$ciiSe[i] < 0) roc$ciiSe[i] <- 0
    if (roc$ciiSe[i] > 1) roc$ciiSe[i] <- 1
    if (roc$cisSe[i] < 0) roc$cisSe[i] <- 0
    if (roc$cisSe[i] > 1) roc$cisSe[i] <- 1
    if (roc$ciiSp[i] < 0) roc$ciiSp[i] <- 0
    if (roc$ciiSp[i] > 1) roc$ciiSp[i] <- 1
    if (roc$cisSp[i] < 0) roc$cisSp[i] <- 0
    if (roc$cisSp[i] > 1) roc$cisSp[i] <- 1
  }

  cat("Numero de individuos medidos:", enf+san, "\n")
  return(roc)
}

set.seed(285899)
IC_SeSp(calculo_SeSp(casos = rnorm(50, 7, 2), controles = rnorm(50, 4, 2)), 50, 50, 0.05)



###############################################################################
# Funcion auxiliar para hacer la gráfica de la curva ROC

plot_roc <- function(roc, plot.ic, auc){
  
  #############################################################################
  # roc: objeto creado por la funcion calculo_SeSp o similares. Contiene datos 
  #      de la sensibilidad y especificidad
  # plot.ic: booleano. Indica si se quiere graficar el IC
  #############################################################################
  
  if(!is.data.frame(roc))   stop("Argumento roc no valido")
  if(!is.logical(plot.ic))  stop("plot.ic debe ser un booleano")
  
  roc <- roc[nrow(roc):1,]
  plot(1-roc$Sp, roc$Se, 
       xlab = "1-Especificidad", ylab = "Sensibilidad", main = "Curva ROC",
       xlim = c(0, 1), ylim = c(0, 1), type = "l")
  abline(a=0, b=1, col= "grey")
  
  # Si pide graficar el intervalo
  if(plot.ic == TRUE){
    lines(1-roc$ciiSp, roc$ciiSe, col = "blue")
    lines(1-roc$cisSp, roc$cisSe, col = "blue")
  }
  
  # texo con el valor del AUC si se ha pedido
  if(is.numeric(auc)){
    texto <- paste("AUC = ", round(auc[1],3))
    text(0.60,0.2,texto, adj = c(0,0))
    if(length(auc)==3){
      texto <- paste("IC = (",round(auc[2],3),",",round(auc[3],3),")")
      text(0.60,0.05,texto, col = "blue", adj = c(0,0))
    }
  }
}

set.seed(285899)
plot_roc(calculo_SeSp(casos = rnorm(50, 7, 2), controles = rnorm(50, 4, 2)), FALSE, auc = FALSE)
plot_roc(IC_SeSp(calculo_SeSp(casos = rnorm(50, 7, 2), controles = rnorm(50, 4, 2)), 50,50, 0.05), TRUE, auc = c(0.864, 0.69598,0.86703))
plot_roc(IC_SeSp(calculo_SeSp(casos = rnorm(500, 7, 1.5), controles = rnorm(300, 5, 2)), 500,300, 0.05), TRUE, auc = 0.77769)



###############################################################################
# Función para el cálculo de curvas ROC no parametricas 

Roc <- function(casos, controles, smooth = FALSE, IC = FALSE, 
                plot = TRUE, plot.ic = FALSE, alpha = 0.05,
                AUC = TRUE, AUC.ic = FALSE, Nsim = 500, pAUC = FALSE,
                descriptivo = FALSE){
  
  #############################################################################
  # casos: Valores de los individuos enfermos
  # Controles: Valores de los individuos sanos
  # smooth: guarda los datos suavizados en vez de los empiricos
  # IC: calcula los intervalos de confianza para la Se y la Sp
  # plot: si plot=TRUE realiza el gráfico de la curva
  # plot.ic: gráfico con los intervalos de confianza o no
  # alpha: nivel de confianza de las estimaciones
  # AUC: si AUC=TRUE calcula el area bajo la curva 
  # AIC.ic: calculo del intervalo de confianza del AUC mediante bootstrap
  # Nsim: numero de simulaciones de bootstrap
  # Descriptivo: grafico descriptivo de las poblaciones de sanos y enfermos
  #############################################################################
  
  # Comprobación de los argumentos
  if(!is.logical(smooth))    stop("El argumento smooth debe ser TRUE o FALSE")
  if(!is.logical(IC))        stop("El argumento IC debe ser TRUE o FALSE")
  if(!is.logical(plot))      stop("El argumento plot debe ser TRUE o FALSE")
  if(!is.logical(plot.ic))   stop("El argumento plot.ic debe ser TRUE o FALSE")
  
  #Eliminar los missing
  if(sum(is.na(casos))!=0){
    warning("Los missing de los casos se eliminarán")
    casos <- casos[!is.na(casos)] }
  if(sum(is.na(controles))!=0){
    warning("Los missing de los controles se eliminarán")
    controles <- controles[!is.na(controles)] }
  # OTros warnings
  if(plot == FALSE & plot.ic == TRUE){
    warning("No se ha soliciado el gráfico principal.Poner plot = TRUE")}
  
  
  # Descriptivo de las poblaciones
  if (descriptivo == 1) poblaciones1(casos, controles)
  if (descriptivo == 2) poblaciones2(casos, controles)  
  
  # Obtencion de la sensibilidad y la especificidad
  enf <- length(casos)
  san <- length(controles)
  n <- enf + san
  # suavizado o no de los datos
  if(smooth == TRUE){
    roc <- calculo_SeSp_suavizado(casos, controles)
  }else{ 
    roc <- calculo_SeSp(casos, controles) }
  
  # Calculo del IC de la curva ROC para la curva suponiendo normalidad
  if(IC == TRUE){
    roc <- IC_SeSp(roc, enf, san, alpha)
  }

  # Cálculo del AUC
  auc <- NULL
  if(AUC[1]==TRUE & length(AUC)==1){
    auc <- calculo_auc(roc)
    if(AUC.ic == TRUE) {
      auc <- c(auc, ICauc_bootstrap(casos, controles, Nsim, alpha)) # PFFF
    }
  }
  
  # Calculo del pauc
  pauc <- NULL  
  if(is.numeric(pAUC) & length(pAUC)==2){
    pauc <- calculo_pauc(roc, pAUC)
  }
  
  # Dibujo de la curva ROC
  if(plot==TRUE) plot_roc(roc, plot.ic, auc)
  
    return(list(casos, controles, roc, auc, pauc)) # auc
}


# Roc(casos, controles, smooth = FALSE, IC = FALSE, plot = TRUE, plot.ic = FALSE,
# alpha = 0.05, AUC = TRUE, AUC.ic = FALSE, Nsim = 300, descriptivo = FALSE)
set.seed(285899)
Roc(casos = rnorm(50, 7, 1.5), controles = rnorm(50, 4, 2)) # 0.9085
roc <- Roc(casos = rnorm(50, 7, 1.5), controles = rnorm(50, 4, 2), plot.ic = TRUE, IC = TRUE, AUC.ic = TRUE, descriptivo = 2)
Roc(casos = rnorm(50, 7, 1.5), controles = rnorm(50, 4, 2),
    plot.ic = TRUE, IC = TRUE, AUC.ic = TRUE, smooth = TRUE)
Roc(casos = rnorm(50, 7, 1.5), controles = rnorm(50, 4, 2), pAUC = c(0.6,0.8))

# library(pROC)
# ?roc
# plot(roc(cases=casos, controls=controles))
# roc(cases=0:10, controls= -5:5)



###############################################################################
# Funcion auxiliar para el calculo de la sensibilidad y especificidad tomando 
# una distribucion conocida para sanos y enfermos.

calculo_SeSp_parametrico <- function(casos, controles, puntos){
  
  #############################################################################
  # casos: Valores de los parametros de la distribucion de enfermos
  # Controles: Valores de los parametros de la distribucion de sanos
  # puntos: numero de divisiones a realizar
  #############################################################################
  
  if(!is.numeric(casos) | length(casos)!= 2){
    stop("Hay que porporcionar los parámetros de la distribución de los casos") }
  if(!is.numeric(controles)| length(controles)!= 2){
    stop("Hay que porporcionar los parámetros de la distribución de los controles") }
  if(!is.numeric(puntos)) stop("El argumento puntos debe toomar un valor numérico")
  
  cat(" Se asumen distribuciones normales para casos y controles. ", "\n")
  # media-5*sd = 99.99 probabilidad
  # minimo y maximo para el dominio de valores factibles
  con_min <- casos[1]-5*casos[2] 
  cas_max <- controles[1]+5*controles[2]
  cortes <- c(-Inf, seq(from = con_min, to = cas_max, by = (cas_max - con_min)/puntos), Inf)
  roc <- data.frame(cortes)
  
  # Cálculo de sensibilidad y especificidad para diferentes puntos de corte
  roc$Se[1] <- 1
  roc$Sp[1] <- 0
  roc$Se[length(cortes)] <- 0
  roc$Sp[length(cortes)] <- 1 
  for(i in 2:length(cortes)-1){
    # repasar esto
    roc$Se[i] <- pnorm(cortes[i], casos[1], casos[2], lower.tail = FALSE)
    roc$Sp[i] <- pnorm(cortes[i], controles[1], controles[2], lower.tail = TRUE)
  }
  
  return(roc)
}

set.seed(1996)
calculo_SeSp_parametrico(casos = c(3,1), controles = c(1,1), puntos=100)



###############################################################################
# Funcion auxiliar para hacer un gráfico de las distintas poblaciones

poblaciones3 <- function(casos, controles){
  
  #############################################################################
  # casos: Valores de los individuos enfermos
  # Controles: Valores de los individuos sanos
  #############################################################################
  
  if(!is.numeric(casos))     stop("Los casos deben ser numericos")
  if(!is.numeric(controles)) stop("Los controles deben ser numericos") 
  
  curve(dnorm(x,casos[1],casos[2]),col = "red",
        xlim = c(min(controles[1]-5*controles[2],casos[1]-5*casos[2]), max(casos[1]+5*casos[2],controles[1]+5*controles[2])),
        ylim = c(0,0.4/min(casos[2], controles[2])),
        main = "Poblaciones", xlab = "Valor del biomarcador", ylab = "")
  curve(dnorm(x,controles[1],controles[2]), add=TRUE)
}

set.seed(285899)
poblaciones3( casos = c(7, 1), controles = c(4, 2))



###############################################################################
# Función para el cálculo de curvas ROC parametricas teniendo la distribucion 
# de las poblaciones de sanos y enfermos. Se momeno la distribucion de los 
# sanos y enfermos solo puede ser normal.

Roc_parametrica <- function(casos, controles, 
                            AUC = TRUE, plot = TRUE,  alpha = 0.05, 
                            descriptivo = FALSE ,puntos= 516){
  
  #############################################################################
  # casos : parametros de la distribucion de los enfermos
  # controles : parametros de la distribucion de los sanos
  # dist.casos : distribucion de los enfermos 
  # dist.controles : distribucion de los sanos
  # plot: si plot=TRUE realiza el gráfico de la curva
  # AUC: si AUC=TRUE calcula el area bajo la curva 
  # alpha: nivel de confianza para las estimaciones
  # puntos: número de divisiones del espacio para el calculo de la aproximacion 
  # descriptivo: indica si hay que graficar las poblaciones de sanos y enfermos   
  #############################################################################
  
  # Comprobación de los argumentos
  if(!is.numeric(casos) | length(casos)>4){
    stop("Hay que porporcionar los parámetros de la distribución de los casos") }
  if(!is.numeric(controles)| length(controles)>4){
    stop("Hay que porporcionar los parámetros de la distribución de los controles") }
  
  roc <- calculo_SeSp_parametrico( casos, controles, puntos)
  
  # Descriptivo de las poblaciones
  if (descriptivo == TRUE) poblaciones3(casos, controles)
  
  # Cálculo del AUC
  auc <- NULL
  if(AUC[1]==TRUE & length(AUC)==1){
    auc <- calculo_auc(roc)
  }
  
  # Calculo del pauc
  pauc <- NULL  
  if(is.numeric(AUC) & length(AUC)==2){
    pauc <- calculo_pauc(roc, AUC)
  }
  
  # Dibujo de la curva ROC
  if(plot==TRUE)  plot_roc(roc, plot.ic = FALSE, auc)
  
  return(list(casos, controles, roc, auc, pauc))
}

Roc_parametrica(casos =  c(185,40), controles =  c(60,50), descriptivo = TRUE)




###############################################################################
# Funcion auxiliar para  el calculo de funcion PruebasBinarias() en la funcion
# PCO(). Solo se puede aplicar en el estudios de casos y controles, por tanto
# hay que introducir la prevalencia

PB <- function(roc, pco, prev, alpha, n){
  
  #############################################################################
  # roc: Objeto creado por la función Roc o Roc_parametrica
  # pco: punto de corte para el que se calcula las medidas
  # prev: nivel de prevalencia
  # alpha: nivel de confianza de las estimaciones
  # n: variable auxilar para el calculo con Roc_parametricas
  #############################################################################
  
  dico <- list()
  for (i in 1:length(pco)){
    # Para las roc
    if(length(roc[[1]])>2){
      casos <- roc[[1]]
      controles <- roc[[2]]
      VP <- as.numeric(sum(casos >= pco[i]))
      VN <- as.numeric(sum(controles < pco[i]))
      FP <- as.numeric(sum(casos < pco[i]))
      FN <- as.numeric(sum(controles >= pco[i]))
      dico[[i]] <- PruebasBinarias(VP, VN, FP, FN, alpha, prev, FALSE)
    }
    
    # para roc_parametrica
    if(length(roc[[1]])==2){
      casos <- roc[[1]]
      controles <- roc[[2]]
      VP <- round(pnorm(pco[i], casos[1], casos[2], lower.tail = FALSE)*n)
      VN <- round(pnorm(pco[i], controles[1], controles[2], lower.tail = TRUE)*n)
      FP <- round(pnorm(pco[i], controles[1], controles[2],lower.tail = FALSE)*n)
      FN <- round(pnorm(pco[i], casos[1], casos[2], lower.tail = TRUE)*n)
      dico[[i]] <- PruebasBinarias(VP, VN, FP, FN, alpha, prev, FALSE )
    }
  }
  return(dico)
}

PB(Roc(rnorm(50,6,1), rnorm(50,3,1)), pco = 4.5, prev = 0.2, alpha = 0.05, n = 10000)
PB(Roc_parametrica(c(6,1), c(3,1)), pco = c(4.5, 5), prev = 0.2, alpha = 0.05, n = 10000)



###############################################################################
# Función para elecció de un punto de corte según un criterio de interés. Hay
# 4 criterios diferentes implementados.

PCO <- function( roc, metodo, prev = NULL, alpha = 0.05, casos = 10000){
  
  #############################################################################
  # roc: data.frame creado por la funcion Roc(). Punto de corte junto con la 
  #      sensibilidad y especificidad asociados
  # metodo: criterio de seleccion del punto de corte óptimo
  # prev: Indica si se quiere calcular la prueba binaria. por defecto es null
  #      y no lo calcula. Si se quieren calcular hay que introducir el valor 
  #      de la prevalencia.
  # alpha: nivel de confianza para las estimaciones.
  # casos: variable auxiliar para calcular las estimaciones en casos de curvas
  #      roc parametricas. No se tiene el numero de casos, por lo que se 
  #      hace proporcional a la distribucion . Si se quiere controlar la 
  #      semilla hay que hacerlo fuera del programa.
  #############################################################################
  
  # if(!is.list(roc)) stop("Argumento roc no valido")
  if(metodo!="youden" &  metodo!="(0,1)" & metodo!="maxSe" & metodo!="maxSp"){ 
    stop("Metodo no valido.")}
  if(alpha<0 | alpha>1)   stop("alpha tiene que ser mayor que 0 y menor que 1")
  # if(is.null(prev) | (prev<0 | prev>1)) stop("prev tiene que ser mayor que 0 y menor que 1")
  
  aux <- roc
  roc <- roc[[3]]
  
  # Youden
  if(metodo == "youden"){
    roc$youden <- roc$Se - (1 - roc$Sp)
    pco <- roc$cc[roc$youden == max(roc$youden)] 
    roc <- roc[ ,!colnames(roc)=="youden"]
  }
  
  # Punto mas cercano a (0,1)
  if(metodo == "(0,1)"){
    roc$dist_cerouno <- sqrt((1-roc$Sp)^2+(roc$Se-1)^2)
    pco <- roc$cc[roc$dist_cerouno == min(roc$dist_cerouno)]
    roc <- roc[ ,!colnames(roc)=="dist_cerouno"]
  }
  
  # maximizar sensibilidad
  if(metodo == "maxSe"){
    pco <- max(roc$cc[roc$Se == 1])
  }
  
  # maximizar la especifidad
  if(metodo == "maxSp"){
    pco <- min(roc$cc[roc$Sp == 1])
  }
  
  PB <- NULL
  if(!is.null(prev)){
    PB <- PB(aux, pco, prev, alpha, casos)
  }
  
  if(!is.null(prev)){
    return(list(pco, PB))
  }else{
    return(pco)
  }
}
  

set.seed(400)
p1 <- Roc(casos = rnorm(50, 7, 2), controles = rnorm(50, 4, 2))
PCO(p1, metodo = "youden") # 4.962894 5.634625 5.796323
PCO(p1, metodo = "(0,1)")  # 4.962894 5.634625
PCO(p1, metodo = "maxSe")  # 1.975445
PCO(p1, metodo = "maxSp")  # 9.536409

PCO(p1, metodo = "(0,1)", prev = 0.1)

# library(OptimalCutpoints)
# ?optimal.cutpoints
# data <- data.frame(x = c(casos, controles) ,status = c(rep(1,50),rep(0,50)), id = 1:100)
# cc1 <- optimal.cutpoints("x" ,"status" , tag.healthy = 0, methods = "Youden",  data)
# cc2 <- optimal.cutpoints("x" ,"status" , tag.healthy = 0, methods = "ROC01" ,  data)
# cc3 <- optimal.cutpoints("x" ,"status" , tag.healthy = 0, methods = "MaxSe",  data)
# cc4 <- optimal.cutpoints("x" ,"status" , tag.healthy = 0, methods = "MaxSp",  data)


