###############################################################################
# Funciones probadas
###############################################################################

###############################################################################
# Cálculo del intervalo de confianza mediante una distribución binomial dado
# el tamaño de la muestra y el número de elementos que presentan la 
# característica de interés. Función auxiliar para la estimación binomial para
# la sensibilidad, la especificidad tanto en casos controles como en estudios
# de cohortes y para el VPP y VPN en estudios de cohortes.

ICbin <- function(x, m, alpha){
  
  #############################################################################
  # x : número de elementos de la muestra con la característida de interés.
  # m : número de elementos de la muestra.
  # alpha : nivel de confianza de las estimaciones.
  #############################################################################

  p <- x/m
  inf <- (1+(m-x+1)/(x*qf(1-alpha/2, 2*x, 2*(m-x+1), lower.tail = FALSE)))^(-1)
  sup <- (1+(m-x)/((x+1)*qf(alpha/2, 2*(x+1),2*(m-x), lower.tail = FALSE)))^(-1)
  return(c(p, inf, sup))
}

# ICbin(x, m, alpha)
ICbin(x = 5, m = 10, alpha = 0.05)
ICbin(x = 50, m = 100, alpha = 0.10)


###############################################################################
# Cálculo del intervalo de confianza mediante una distribución normal dado
# el tamaño de la muestra y el número de elementos que presentan la 
# característica de interés. Función auxiliar para la estimación normal para
# la sensibilidad, la especificidad tanto en casos controles como en estudios
# de cohortes y para el VPP y VPN en estudios de cohortes.

ICnorm <- function(a, t, alpha){
  
  #############################################################################
  # a : número de elementos de la muestra con la característida de interés.
  # t : número de elementos de la muestra.
  # alpha : nivel de confianza de las estimaciones.
  #############################################################################
  
  p <- a/t
  inf <- a/t + qnorm(alpha/2)/t*sqrt(a*(t-a)/t)
  sup <- a/t + qnorm(1-alpha/2)/t*sqrt(a*(t-a)/t)
  result <- c(p, inf, sup)
  return(result)
}

# ICnorm(a, t, alpha)
ICnorm(a = 5, t = 10, alpha = 0.10)
ICnorm(a = 50, t = 100, alpha = 0.05)


###############################################################################
# Función auxiliar para el cálculo del intervalo de confianza para el DLR 
# positivo mediante la aproximación a una normal.

ICDLRpnorm <- function(Se, Sp, enf, san, alpha){
  
  #############################################################################
  # Se : Sensibilidad de la prueba.
  # Sp : Especificidad de la prueba.
  # enf : número de elementos de la muestra con la caracteristica de interés.
  # san : número de elementos de la muestra sin la caracteristica de interés.
  # alpha : nivel de confianza de las estimaciones.
  #############################################################################
  
  p <- Se/(1-Sp)
  inf <- exp(log(Se/(1-Sp))+qnorm(alpha/2)*sqrt((1-Se)/enf/Se)+Sp/san/(1-Sp))
  sup <- exp(log(Se/(1-Sp))+qnorm(1-alpha/2)*sqrt((1-Se)/enf/Se)+Sp/san/(1-Sp))
  result <- c(p, inf, sup)
  return(result)
}

# ICDLRpnorm(Se, Sp, enf, san, alpha)
ICDLRpnorm(Se = 0.9, Sp = 0.8, enf = 40, san = 80, alpha = 0.05)
ICDLRpnorm(Se = 0.8, Sp = 0.9, enf = 50, san = 80, alpha = 0.10)
    

###############################################################################
# Función auxiliar para el cálculo del intervalo de confianza para el DLR 
# negativo mediante la aproximación a una normal.

ICDLRnnorm <- function(Se, Sp, enf, san, alpha){
  
  #############################################################################
  # Se : Sensibilidad de la prueba.
  # Sp : Especificidad de la prueba.
  # enf : número de elementos de la muestra con la caracteristica de interés.
  # san : número de elementos de la muestra sin la caracteristica de interés.
  # alpha : nivel de confianza de las estimaciones.
  #############################################################################
  
  p <- (1-Se)/Sp
  inf <- exp(log((1-Se)/Sp)+qnorm(alpha/2)*sqrt((1-Se)/enf/Se)+Sp/san/(1-Sp))
  sup <- exp(log((1-Se)/Sp)+qnorm(1-alpha/2)*sqrt((1-Se)/enf/Se)+Sp/san/(1-Sp))
  result <- c(p,inf,sup)
  return(result)
}

# ICDLRnnorm(Se, Sp, enf, san, alpha)
ICDLRnnorm(Se = 0.9, Sp = 0.8, enf = 50, san = 50, alpha = 0.05)
ICDLRnnorm(Se = 0.8, Sp = 0.9, enf = 50, san = 50, alpha = 0.10)



###############################################################################
# Función auxiliar para el cálculo del intervalo de confianza para el VPP en 
# estudios de casos  controles

ICVPP_cc <- function(prev, Se, Sp, enf, san, alpha){
  
  #############################################################################
  # prev: prevalencia del evento de interés
  # Se : Sensibilidad de la prueba.
  # Sp : Especificidad de la prueba.
  # enf : número de elementos de la muestra con la caracteristica de interés.
  # san : número de elementos de la muestra sin la caracteristica de interés.
  # alpha : nivel de confianza de las estimaciones.
  #############################################################################
  
  p <- prev*Se/(prev*Se+(1-prev)*(1-Sp))
  DLRp <-ICDLRpnorm(Se, Sp, enf, san, alpha)
  inf <- DLRp[2]*(prev/(1-prev))/(1+(DLRp[2]*(prev/(1-prev))))
  sup <- DLRp[3]*(prev/(1-prev))/(1+(DLRp[3]*(prev/(1-prev))))
  result <- c(p, inf, sup)
  return(result)
}

# ICVPP_cc(prev, Se, Sp, enf, san, alpha)
ICVPP_cc(prev = 0.50, Se = 0.8, Sp = 0.9, enf = 25, san = 25, alpha = 0.05)
# OJO, intervalo no centrado si no tiene congruencia con los datos
ICVPP_cc(prev = 0.10, Se = 0.9, Sp = 0.9, enf = 100, san = 100, alpha = 0.05)



###############################################################################
# Función auxiliar para el cálculo del intervalo de confianza para el VPN en 
# estudios de casos  controles

ICVPN_cc <- function(prev, Se, Sp, enf, san, alpha){
  
  #############################################################################
  # prev: prevalencia del evento de interés
  # Se : Sensibilidad de la prueba.
  # Sp : Especificidad de la prueba.
  # enf : número de elementos de la muestra con la caracteristica de interés.
  # san : número de elementos de la muestra sin la caracteristica de interés.
  # alpha : nivel de confianza de las estimaciones.
  #############################################################################
  
  p <- (1-prev)*Sp/((1-prev)*Sp+prev*(1-Se))
  DLRn <-ICDLRnnorm(Se, Sp, enf, san, alpha)
  inf <- DLRn[2]*(prev/(1-prev))/(1+(DLRn[2]*(prev/(1-prev))))
  sup <- DLRn[3]*(prev/(1-prev))/(1+(DLRn[3]*(prev/(1-prev))))
  result <- c(p, inf, sup)
  return(result)
}

# ICVPN_cc(prev, Se, Sp, enf, san, alpha)
ICVPN_cc(prev = 0.50, Se = 0.8, Sp = 0.9, enf = 25, san = 25, alpha = 0.05)
ICDLRnnorm(Se = 0.8, Sp = 0.9, enf = 25, san = 25, alpha = 0.05)
# OJO, intervalo no centrado  si no tiene congruencia con los datos
ICVPN_cc(prev = 0.10, Se = 0.9, Sp = 0.9, enf = 100, san = 100, alpha = 0.05)


###############################################################################
# Función que calcula la estimacion de la sensibilidad, especificidad, VPP, 
# VPN, DLR+ y DLR- de una prueba diagnóstica con resultado binario. Puede 
# calcular mediante una aproximacion normal o la distribución binomial o tener
# en cuenta si el estudio es de cohortes o de casos y controles.

PruebasBinarias <- function(VP, VN, FP, FN, alpha=0.05, prev=-1, aprox=FALSE){
  
  #############################################################################
  # VP : número de verdaderos positivos
  # VN : número de verdaderos negativos
  # FP : número de falsos positivos
  # FN : número de falsos negativos
  # alpha : nivel de confianza de las estimaciones
  # prev: prevalencia. Solo se necesita para estimar en casos y controles
  #############################################################################
  
  if(alpha <= 0 | alpha >= 1){
    stop("Argumento alpha no válido, debe ser mayor que cero e inferior a 1") }
  if(VP <= 0 | VN <= 0 | FP <= 0 | FN <= 0){
    stop("VP, VN, FP, FN son frecuencias, deben ser mayores a 0")}
  if(is.integer(VP) | is.integer(VN) | is.integer(FP) | is.integer(FN)){
    stop("VP, VN, FP, FN deben ser valores enteros") }
  if(prev != -1 & (prev <= 0 | prev >= 1)){
    stop("La prevalencia debe ser > 0 e < 1") }
  
  # Variables auxiliares:
  n   <- VP+VN+FP+FN # número total de individuos estudiados
  enf <- VP+FN       # número total de enfermos
  san <- n-enf       # número total de sanos
  pos <- VP+FP       # número total de resultados positivos
  neg <- n-pos       # número total de resultados negativos
  
  # Estimaciones binomiales:
  if( aprox == FALSE ){
    Se <- ICbin(VP, enf, alpha)
    Sp <- ICbin(VN, san, alpha)
    if (prev == -1){
      # Estudio de cohortes
      VPP <- ICbin(VP, pos, alpha)
      VPN <- ICbin(VN, neg, alpha)
    }else{
      # Estudio de Casos-Control
      VPP <- ICVPP_cc(prev, Se[1], Sp[1], enf, san, alpha)
      VPN <- ICVPN_cc(prev, Se[1], Sp[1], enf, san, alpha) }
    
  }else{
    
    # Estimaciones normales:
    Se <- ICnorm(VP, enf, alpha)
    Sp <- ICnorm(VN, san, alpha)
    if (prev == -1){
      # Estudio de cohortes
      VPP <- ICnorm(VP, pos, alpha)
      VPN <- ICnorm(VN, neg, alpha)
    }else{
      # Estudio de Casos-Control
      VPP <- ICVPP_cc(prev, Se[1], Sp[1], enf, san, alpha)
      VPN <- ICVPN_cc(prev, Se[1], Sp[1], enf, san, alpha) }
  }
  # Para el DLR solo hay estimación normal
  DLRp <- ICDLRpnorm(Se[1], Sp[1], enf, san, alpha)
  DLRn <- ICDLRnnorm(Se[1], Sp[1], enf, san, alpha)
  result <- as.data.frame(rbind(Se, Sp, VPP, VPN, DLRp, DLRn))
  colnames(result) <- c("p","inf","sup")
  
  cat("---------------------------------------------------------","\n")
  if(prev==-1) {cat("Estudio de cohortes","\n")
  }else{cat("Estudio de casos-control","\n")}
  if(prev==-1) prev <- enf/n
  cat("Prevalencia:",prev*100,"%","\n")
  cat("Nivel de confianza:", 100-(100*alpha),"%","\n")
  cat("---------------------------------------------------------","\n")
  cat("                   Prueba de referencia","\n")
  # Falta decir cual es la categoría de la prueba de referencia o TOTAL.
  data <- data.frame(matrix(c("1", "0", "Total", VP, FN, enf, FP, VN, san, pos, neg, n), 3))
  colnames(data) <- c("      ","     1","     0"," Total")
  rownames(data) <- c("  Prueba", "diagnostica","")
  print(data)
  cat("---------------------------------------------------------","\n")
  
  return(round(result,6))
}

# PruebasBinarias(VP, VN, FP, FN, alpha=0.05, prev=-1, aprox=FALSE)
# PruebasBinarias(164,85,0,33)
PruebasBinarias(164,85,18,33)
PruebasBinarias(164,85,18,33, aprox = TRUE)
PruebasBinarias(164,85,18,33, prev = 0.1) # VPN
PruebasBinarias(164,85,18,33, aprox = TRUE, prev = 0.1)
PruebasBinarias(50,50,50,50)
PruebasBinarias(100,100,1,1) # OJO CON LOS DLR



### Pruebas binarias

###############################################################################
# Función que calcula la estimacion de la sensibilidad, especificidad, 
# de una prueba diagnóstica con resultado binario comparando con una prueba de 
# referencia imperfecta. Se puede calcular mediante una aproximacion normal o
# la distribución binomial o tener en cuenta si el estudio es de cohortes o de
# casos y controles.

## Prueba de referencia imperfecta:
PruebasBinariasImperfectas <- function(VP, VN, FP, FN, Se1, Sp1, alpha = 0.05, aprox=FALSE, prev = -1){
  
  #############################################################################
  # VP : número de verdaderos positivos
  # VN : número de verdaderos negativos
  # FP : número de falsos positivos
  # FN : número de falsos negativos
  # Se1 : Sensibilidad de la prueba de referencia
  # Sp1 : Especificidad de la prueba de referencia
  # alpha : nivel de confianza de las estimaciones
  #############################################################################
  
  if(alpha <= 0 | alpha >= 1){
    stop("Argumento alpha no válido, debe ser mayor que cero e inferior a 1") }
  if(Se1 <= 0 | Se1 >= 1 ){
    stop("Argumento Se1 no válido, debe ser mayor que cero e inferior a 1") }
  if(Sp1 <= 0 | Sp1 >= 1 ){
    stop("Argumento Sp1 no válido, debe ser mayor que cero e inferior a 1") }
  if(VP <= 0 | VN <= 0 | FP <= 0 | FN <= 0){
    stop("VP, VN, FP, FN son frecuencias, deben ser mayores a 0")}
  if(is.integer(VP) | is.integer(VN) | is.integer(FP) | is.integer(FN)){
    stop("VP, VN, FP, FN deben ser valores enteros") }
  if(prev != -1 & (prev <= 0 | prev >= 1)){
    stop("La prevalencia debe ser > 0 e < 1") }
  
  # Variables auxiliares:
  n   <- VP+VN+FP+FN # número total de individuos estudiados
  enf <- VP+FN       # número total de enfermos
  san <- n-enf       # número total de sanos
  pos <- VP+FP       # número total de resultados positivos
  neg <- n-pos       # número total de resultados negativos
  
  # Estimaciones binomiales:
  if( aprox == FALSE ){
    Se <- ICbin(round(((VP + FP) * Sp1 - FP)/((VP + FN) - (1 - Sp1) * n)*enf,4), enf, alpha)
    Sp <- ICbin(round(((FN + VN) * Se1 - FN)/(n * Se1 - VP - FN)*san, 4), san, alpha)
    
  }else{
    # Estimaciones normales:
    Se <- ICnorm(round(((VP + FP) * Sp1 - FP)/((VP + FN) - (1 - Sp1) * n)*enf, 4), enf, alpha)
    Sp <- ICnorm(round(((FN + VN) * Se1 - FN)/(n * Se1 - VP - FN)*san, 4), san, alpha)
  }
  
  # Para el DLR solo hay estimación normal
  DLRp <- ICDLRpnorm(Se[1], Sp[1], enf, san, alpha)
  DLRn <- ICDLRnnorm(Se[1], Sp[1], enf, san, alpha)
  
  result <- as.data.frame(rbind(Se, Sp, DLRp, DLRn)) #, VPP, VPN
  colnames(result) <- c("p","inf","sup")
  
  cat("---------------------------------------------------------","\n")
  if(prev==-1) prev <- (((VP + FP) / n + Sp1 - 1)) / ( Se1 + Sp1 - 1)
  if(prev==-1) {cat("Estudio de cohortes","\n")
  }else{cat("Estudio de casos-control","\n")}
  cat("Prevalencia:",prev*100,"%","\n")
  cat("Nivel de confianza:", 100-(100*alpha),"%","\n")
  cat("Sensibilidad de la prueba de referencia :", Se1 ,"\n")
  cat("Especificidad de la prueba de referencia :", Sp1 ,"\n")
  cat("---------------------------------------------------------","\n")
  cat("            Prueba de referencia imperfecta","\n")
  data <- data.frame(matrix(c("1", "0", "Total", VP, FN, enf, FP, VN, san, pos, neg, n), 3))
  colnames(data) <- c("      ","     1","     0"," Total")
  rownames(data) <- c("  Prueba", "diagnostica","")
  print(data)
  cat("---------------------------------------------------------","\n")
  
  return(result)
}
PruebasBinariasImperfectas(164,85,18,33, Se1 = 0.95, Sp1 = 0.80)



###############################################################################
# Función para el cálculo del tamaño muestral necesario para alcanzar la
# pontencia y el nivel de confianza deseado. Unicamente se calcula para
# pruebas binarias en estudios de casos y controles.

TamanoMuestral <- function(FPF1, TPF1, TPF0, FPF0, alpha = 0.05, beta = 0.2){
  
  #############################################################################
  # FPF1 : 1-Especificidad deseada
  # TPF1 : Sensibilidad deseada
  # FPF0 : 1-Especificidad mínima aceptable
  # TPF0 : Sensibilidad mínima aceptable
  # alpha :  nivel de confianza del contraste del que se obtiene el tamaño 
  #        muestral
  # beta : 1 - potencia del contraste del que se obtiene el tamaño muestral
  #############################################################################
  
  nE <- ((qnorm(sqrt(1-alpha))*sqrt(TPF0*(1-TPF0)))+(qnorm(sqrt(1-beta))*sqrt(TPF1*(1-TPF1))))**2/(TPF1-TPF0)**2
  nS <- ((qnorm(sqrt(1-alpha))*sqrt(FPF0*(1-FPF0)))+(qnorm(sqrt(1-beta))*sqrt(FPF1*(1-FPF1))))**2/(FPF1-FPF0)**2

  cat( "\n")   
  cat("Para alcanzar una potencia de",(1-beta)*100,"% con un nivel ", "\n")
  cat("de confianza de",(1-alpha)*100,"% se necesita una muestra de", "\n")
  cat( ceiling(nS),"sanos y",ceiling(nE),"enfermos.", "\n") 
  cat( "\n") 
  return(list(ceiling(nE),ceiling(nS)))
}

TamanoMuestral( FPF1 = 0.05, TPF1 = 0.90, TPF0 = 0.75,
                FPF0 = 0.2, alpha = 0.1, beta = 0.1)
