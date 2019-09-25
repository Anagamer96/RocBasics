# RocBasics
En la librería se recogen funciones para trabajar con pruebas diagnosticas


PruebasBinarias. Esta función devuelve una matriz con los valores puntuales y por intervalo de cada una de las medidas descritas junto con unas salidas por consola para comprobar que los datos se han introducido correctamente. Para algunos de los argumentos se describe posteriormente la base teórica. Los argumentos de la función son:
- “VP”: número de verdaderos positivos que se ha obtenido en la prueba. Haciendo referencia a la Tabla 1 es la celda 
- “VN”: número de verdaderos negativos que se ha obtenido en la prueba. Haciendo referencia a la Tabla 1 es la celda         
- “FP”: número de falsos positivos que se ha obtenido en la prueba. Haciendo referencia a la Tabla 1 es la celda        
- “FN”: número de falsos negativos que se ha obtenido en la prueba. Haciendo referencia a la Tabla 1 es la celda 
- “alpha”: nivel de riesgo que se quiere para las estimaciones. Por defecto toma el valor 0.05, lo que corresponde a un nivel de confianza del 95%.
- “prev”: tipo de estudio que se ha realizado, de cohortes o de casos y controles. Por defecto toma el valor -1 lo que corresponde a un estudio de cohortes. Si se quisiera calcular las medidas de eficacia para un estudio de casos y controles habría que introducir el valor estimado de la prevalencia.
- “aprox”: como se ha visto para algunas de las medidas de eficacia hay distintos métodos para el cálculo de los intervalos. Con este argumento se decide si se quieren calcular las estimaciones mediante la distribución exacta o la distribución asintótica. Por defecto toma el valor FALSE lo que indica que se calcula la estimación mediante la distribución exacta.


Para estudiar las medidas de una prueba cuya prueba de referencia es imperfecta se implementa la función PruebasBinariasImperfectas. Los argumentos son:
- “VP”:número de verdaderos positivos que se ha obtenido en la prueba. Haciendo referencia a la Tabla 1 es la celda 
- “VN”: número de verdaderos negativos que se ha obtenido en la prueba. Haciendo referencia a la Tabla 1 es la celda 
- “FP”: número de falsos positivos que se ha obtenido en la prueba. Haciendo referencia a la Tabla 1 es la celda 
- “FN”: número de falsos negativos que se ha obtenido en la prueba. Haciendo referencia a la Tabla 1 es la celda 
- “Se1”: valor de la sensibilidad de la prueba de referencia.
- “Sp1”: valor de la especificidad de la prueba de referencia
- “alpha”: nivel de riesgo que se quiere para las estimaciones. Por defecto toma el valor 0.05, lo que corresponde a un nivel de confianza del 95%.
- “prev”: como se ha dicho previamente la principal diferencia entre los estudios de casos y controles y los estudios de cohortes reside en la estimación de la prevalencia. Por defecto toma el valor -1. Cuando toma este valor la función asume que es un estudio de cohortes. Si por el contrario se quisieran calcular las medidas de eficacia para un estudio de casos y controles se debe introducir el valor estimado de la prevalencia que se tiene de conocimiento a priori.
- “aprox”: indica si se quiere calcular las estimaciones mediante la distribución exacta o asintótica. Por defecto toma el valor FALSE lo que indica que se calcula mediante la distribución exacta.


Roc. Esta función solicita los datos de los enfermos y los sanos que se han recogido en la muestra. Cada uno de estos datos se introduce en los argumentos “casos” y “controles”. Para ampliar la funcionalidad de la función se han incluido diversos argumentos. Para algunos argumentos no se ha explicado la base teórica aun:
-	“casos”: valores del biomarcador obtenido en la muestra de enfermos como un vector.
-	“controles”: valores del biomarcador obtenido en la muestra de sanos como un vector.
-	“smooth”: indica si se quiere realizar un suavizado de los datos para tener una curva en vez de una función escalonada. Por defecto toma el valor FALSE.
-	“IC”: indica si se quiere calcular el intervalo de confianza de la sensibilidad y la especificidad mediante las ecuaciones [13] y [16] para cada punto de corte posible. Por defecto toma el valor FALSE.
-	“plot”: indica si se quiere realizar o no el gráfico de la curva ROC. Por defecto esta codificado como TRUE.
-	“plot.ic” indica si se quiere incluir en el gráfico el intervalo calculado con la función “IC”. Este argumento solo aplica si se ha solicitado el cálculo del intervalo de confianza para la sensibilidad y la especificidad y además se solicita realizar el gráfico. Por defecto esta codificada como FALSE.
-	“AUC”: indica si se quiere calcular el área bajo la curva ROC. Por defecto esta codificado como TRUE. 
-	“pAUC”: indica si se quiere calcular o no el área parcial bajo la curva ROC (pAUC). Por defecto esta codificado como FALSE. Si se quiere calcular se indica mediante un vector de longitud 2 que recoge los valores del complementario de la especificidad entre los que se quiere calcular el área.
-	“AUC.ic”: calcula el intervalo de confianza del AUC utilizando el método de remuestreo de Bootstrap. Por defecto toma el valor FALSE.
-	“Nsim”: número de simulaciones que se quieren realizar en Bootstrap para el AUC.ic. Este parámetro solo aplica si se ha solicitado el cálculo del intervalo de confianza del AUC. Por defecto toma el valor 500.
-	“descriptivo”: indica si se quieren realizar un gráfico descriptivo de las poblaciones de sanos y enfermos muestreados. Por defecto toma el valor FALSE. Hay dos gráficos descriptivos diferentes disponibles. Se indica el valor 1 para que realice dos histogramas en gráficos diferentes situados uno encima del otro o el valor 2 para que realice un único gráfico descriptivo que superponga las distribuciones suavizadas de las poblaciones de sanos y enfermos.


Roc_parametrica. Esta función solo contempla la distribución normal aunque esta no sea la única distribución existente, las demás distribuciones se implementarán en un futuro. Para poner utilizar esta función se necesita conocer los parámetros de las distribuciones de sanos y enfermos, es decir, su media y su desviación típica. Aunque esta función es similar a la función Roc explicada anteriormente algunos argumentos cambian ligeramente mientras que otros desaparecen, por lo que se explican todos de nuevo:
-	“casos”: media y desviación típica de la población de enfermos como un vector de tamaño 2.
-	“controles”: media y desviación típica de la población de sanos como un vector de tamaño 2.
-	“plot”: indica si se quiere realizar o no el gráfico de la curva ROC. Por defecto toma valor TRUE.
-	 “AUC”: indica si se quiere calcular o no el área bajo la curva ROC. Por defecto toma valor TRUE. Si no se quiere calcular hay que indicarlo.
-	“pAUC”: indica si se quiere calcular o no el área parcial bajo la curva ROC (pAUC). Por defecto está codificado como FALSE. Si se quisiera calcular habría que indicarlo mediante un vector de longitud 2 que indicara los valores del complementario de la especificidad entre los que se quiere calcular el área.
-	 “descriptivo”: indica si se quiere realizar o no gráficos descriptivos de las poblaciones de sanos y enfermos muestreados. Si se solicita se realiza un único gráfico con líneas que representan las poblaciones proporcionadas. Este argumento por defecto toma el valor FALSE. 
-	“puntos”: el objetivo de la función es calcular para un rango de valores concreto en muchos puntos equidistantes el valor de la sensibilidad y la especificidad. Este argumento permite seleccionar el número de puntos equidistantes que se quieren emplear. Por defecto está fijado el valor 516 pero se puede seleccionar un valor mayor para tener mayor precisión o menor en caso de querer los resultados más rápido.


PCO halla el punto de corte óptimo y devuelve una lista de dos elementos. En el primer elemento se recogen el punto o los puntos de corte óptimos y en el segundo, si se solicita, se proporcionan las medidas de eficacia asociadas a los puntos de corte obtenidos. Los argumentos son:
-	“roc”: objeto creado por la función Roc o Roc_parametrica.
-	“metodo”: selección del criterio que el usuario desea emplear. Los 4 valores disponibles para este valor son: ‘Youden’, ‘(0,1)’, ‘maxSe’ y ‘maxSp’. Si no se introducen exactamente así los valores (con esas mayúsculas, minúsculas y comas) el argumento no será válido y por tanto no se calculará el punto de corte óptimo.
-	“PB”: indica si se quiere realizar medidas de eficacia de la prueba en caso de que se dicotomizarán los resultados de la prueba mediante el punto de corte óptimo seleccionado. Las medidas se calculan mediante la función PruebasBinarias explicada anteriormente.


TamanoMuestral. Los argumentos de la función son:
-	“FPF1”: especificidad esperada.
-	“TPF1”: sensibilidad esperada.
-	“FPF0”: especificidad mínima aceptable.
-	“TPF0”: sensibilidad mínima aceptable.
-	“alpha”: nivel de riesgo que se asume para el contraste. Por defecto toma el valor 0.05, lo que corresponde a un nivel de confianza del 95%.
-	“beta”: complementario de la potencia del test. Por defecto toma el valor 0.20, lo que corresponde a una potencia de 80%.





