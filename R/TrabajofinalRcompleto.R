 
#1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. Recuerda incluir los nombres de las columnas como el encabezado. Examina el archivo con las funciones head(), tail() y dim(). (1 pto) 

#Tengo que cargar el archivo de texto (que debe estar en la misma #carpeta que este file para poder examinarlo guardo la información del archivo como la variable "datos


datos <- read.table(file.choose(),header = TRUE, sep = "\t", dec = ".")

head(datos) #Muestra las 6 primeras filas
tail(datos) #Muestra las 6 ultimas filas
dim(datos) #Muestra las dimensiones del conjunto de datos (filas y columnas)


### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)


#Sustituimos los puntos de los encabezados por guiones
columna_control <- c("C2.TPM", "C1.TPM", "C3.TPM")
columna_condicion <- c("X1W5.TPM", "X1W6.TPM", "X1W5.COV")


#La media (mean) representa el nivel general de expresion de un gen a lo largo de
#sus 3 réplicas.
#La desviacion estandar (SD)mide la variabiliad entre esas 3 réplicas
#apply(matriz, 1, funcion), el parámetro 1 indica que se debe aplicar
#la funcion recorriendo la tabla fila por fila.

#Calculamos y añadimos las 4 columnas a la tabla "datos"

#Columna 1: Media de réplicas de control
datos$Media_Control <- apply(datos[, columna_control], 1, mean)

#Columna 2: Desviacion estandar de control
datos$SD_Control <- apply(datos[, columna_control], 1, sd)

#Columna 3: Media de réplicas de condición
datos$Media_Condicion <- apply(datos[, columna_condicion], 1, mean)

#Columna 4: Desviación estádar de condición
datos$SD_Condicion <- apply(datos[, columna_condicion], 1, sd)

#Verificacion de que la nueva tabla contiene 11 columnas (11 más 4)
head(datos)



### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

#Asignamos colores a cada variable
color_mediacontrol <- "blue"
color_mediacondicion <-"red"

#utilizamos una funcion para dividir el espacio del gráfico en 1 fila y 2 columnas para
#poder dibujar 2 histogramas dentro de la misma imagen

par(mfrow = c(1, 2))

#Hacemos un histograma que represente visualmente la expresion media de los genes control
#vs condición.

hist(datos$Media_Control,
	col = color_mediacontrol,
	main = "Histograma Media Control",
	xlab = "Nivel de Expresion (TPM)",
	ylab = "Frecuencia de Genes")

hist(datos$Media_Condicion,
	col = color_mediacondicion,
	main = "Histograma Media Condicion",
	xlab = "Nivel de Expresion (TPM)",
	ylab = "Frecuencia de Genes")

#reestablecemos la ventana gráfica a la configuracion normal

par(mfrow = c(1, 1))


### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

#El boxplot representa la mediana, cuartiles y outliers de cada grupo
#Permite comparar facilmente si la condicion provoco un cambio en
#la expresion frente al control.

boxplot(datos$Media_Control, datos$Media_Condicion,
	col = c(color_mediacontrol, color_mediacondicion),
	names = c("Control", "Condicion"),
	main = "Comparcion de Expresión: Control vs Condición",
	ylab = "Valores de Expresión (TPM)")

#añadimos leyenda explicando colores con legend() y la ubicamos en la esquina izq con topleft
legend("topleft",
	legend = c("Control", "Condicion"),
	fill = c(color_mediacontrol, color_mediacondicion))

	


### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

#El test shapiro estudia las medias globales para ver si siguen una distribución normal
#Hacemos el pvalue para determinar si la diferencia entre las replicas del control
#las de la condicion en un gen específico es significativa (real).
shapiro_control <- shapiro.test(datos$Media_Control)$p.value
shapiro_condicion <- shapiro.test(datos$Media_Condicion)$p.value

cat("Shapiro-Wilk p-value (Media Control):", shapiro_control, "/n")
cat("Shapiro-Wilk p-value (Media Condicion):", shapiro_condicion, "/n")

#Aplicamos la funcion gen a gen por dila, guardando el pvalue obtenido
 p_values <- apply(datos, 1, function(fila) {
	ctrl <- as.numeric(fila[columna_control])
	cond <- as.numeric(fila[columna_condicion])
	resultado <- t.test(ctrl, cond)
	return(resultado$p.value) #Este return permite devolver únicamente la cifra pvalue
})

#Añadimos la columna
datos$p_value <- p_values

#Confirmamos que la columna se ha añadido correctamente
head(datos)



### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

#A los genes que tengan un pvalue menor que 0.05 los llamaremos significativos
genes_significativos <- datos[datos$p_value < 0.05, ]

#Contamos el número de genes significativos mediante nrow()
num_significativos <- nrow(genes_significativos)
cat("Numero total de genes significativos expresados (p < 0.05):", num_significativos, "\n")


#De los genes significativos, contamos los sobreexpresados y los reprimidos

if (num_significativos > 0){
	sobreexpresados <- sum(genes_significativos$Media_Condicion > genes_significativos$Media_Control)
	reprimidos <- sum(genes_significativos$Media_Condicion < genes_significativos$Media_Control)


	#Localizamos nombres exactos de los genes significativos
	cat("Genes sobreexpresados (Condición > Control):", sobreexpresados, "\n")
	cat("Genes reprimidos (Condición < Control):", reprimidos, "\n")
	cat("\n--- Listado de genes diferencialmente expresados ---\n")
	print(genes_significativos$Gene.Name)
} else {
	cat("No se obtuvieron genes con p-value <0.005.\n")
}

#Por tanto, tenemos 362 genes significativos (con un p.value < 0.05)
#De los cuales 50 estan sobreexpresados y 312 reprimidos
#Las diferencias observadas entre la expresion de control y condicion de estos genes
#son estadisticamente significativos



### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

#Creamos matriz numerica extrayendo las columnas de las replicas
matriz_heatmap <- as.matrix(genes_significativos[, c(columna_control, columna_condicion)])
	
#Asignamos el nombre de cada gen como leyenda para las filas de la matriz
rownames(matriz_heatmap) <- genes_significativos$Gene.Name

#Generamos un heatmap: scale estandariza los valores por fila (gen) para comparar expresion
#heat.colors(256) aplica colores degradados para visualizar el heatmap
#main añade un título al gráfico
heatmap(matriz_heatmap,
	scale = "row",
	col = heat.colors(256),
	main = "Heatmap de genes significativos")

