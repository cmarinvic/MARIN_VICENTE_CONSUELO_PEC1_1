
# Instalar Bioconductor si no está instalado
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar SummarizedExperiment si no está instalado
BiocManager::install("SummarizedExperiment")

#Instalar pillar
install.packages("pillar")
library(pillar)


# Cargar la librería
library(SummarizedExperiment)



# Leer el archivo de datos de expresión

# Leer el archivo completo como texto
data_matrix <- read.table(
  "C:\\Users\\Usuario\\Desktop\\2020\\UOC\\Master_bioinformatica_bioestadistica\\Analisis_datos_omicos\\PEC1\\ST000567_AN000871_Results_1.txt",
  header = TRUE,    # La primera fila es un encabezado
  sep = "\t",       # Separador de columnas es tabulación
  row.names = 1,    # La primera columna será utilizada como nombres de fila
  
)

data_matrix

str(data_matrix)



# Cargar la librería SummarizedExperiment si no la tienes cargada
library(SummarizedExperiment)

# Crear el objeto SummarizedExperiment
se_object <- SummarizedExperiment(
  assays = list(counts = data_matrix),  # Los datos de expresión
  # No hay metadatos de las columnas (o puedes añadirlos si los tienes)
)

# Verifica el objeto SummarizedExperiment
se_object



#Inspección del objeto:

# Ver el objeto SummarizedExperiment
se_object

# Obtener dimensiones de la matriz de expresión
dim(se_object)

# Ver los nombres de las filas (metabolitos o características)
head(rownames(se_object))

# Ver los nombres de las columnas (muestras)
head(colnames(se_object))



#Estadística descriptiva

# Obtener estadísticas descriptivas de los datos (matriz de expresión)
summary(assays(se_object)$counts)  # Si el assay se llama "counts"

#Mean
# Acceder a la matriz de datos (generalmente 'counts')
data_matrix <- assays(se_object)$counts

# Calcular la media por fila excluyendo la primera columna
row_means <- rowMeans(data_matrix[, 2:ncol(data_matrix)], na.rm = TRUE)
row_sds <- apply(assays(se_object)$counts, 1, sd)

# Ver las primeras medias calculadas
head(row_means)



# Calcular la media y desviación estándar por muestra (columna)
col_means <- colMeans(assays(se_object)$counts)
col_sds <- apply(assays(se_object)$counts, 2, sd)




#Visualización de datos:

# Histograma de las intensidades para todos los metabolitos (filas)
hist(log2(row_means + 1), main = "Histograma de intensidades por metabolito", xlab = "Log2 Intensidad")

# Boxplot de las intensidades para las muestras (columnas)
boxplot(assays(se_object)$counts, main = "Boxplot de las intensidades por muestra", las = 2)

#Escalado y normalización de los datos:

# Log-transformación de los datos de expresión (si no se ha hecho previamente)
log_counts <- log2(assays(se_object)$counts + 1)

# Boxplot después de log-transformar los datos
boxplot(log_counts, main = "Boxplot de las intensidades log-transformadas por muestra", las = 2)

#Correlación entre muestras:

# Matriz de correlación entre las muestras
correlation_matrix <- cor(t(assays(se_object)$counts))

# Visualización de la matriz de correlación como un heatmap
library(pheatmap)
pheatmap(correlation_matrix, main = "Matriz de Correlación entre Muestras")


# Realizar PCA sobre los datos (log-transformados si es necesario)
pca_result <- prcomp(t(log_counts), center = TRUE, scale. = TRUE)

# Visualizar los resultados del PCA
library(ggplot2)
pca_data <- data.frame(pca_result$x)

# Graficar las primeras dos componentes principales
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = colData(se_object)$condition)) +  # Cambia 'condition' por el nombre de la variable de interés
  labs(title = "PCA de las Muestras", x = "PC1", y = "PC2") +
  theme_minimal()















#_____________________________________________________________









# Instalar Bioconductor si no está instalado
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar SummarizedExperiment si no está instalado
BiocManager::install("SummarizedExperiment")


# Cargar la librería
library(SummarizedExperiment)

# Ejemplo de datos (puedes reemplazarlo con tus datos reales)
data_matrix <- matrix(rnorm(100), nrow = 10)  # 10 filas y 10 columnas (ejemplo)

# Metadatos de las filas (por ejemplo, información sobre los genes/metabolitos)
row_metadata <- data.frame(GeneID = paste0("Gene", 1:10))

# Metadatos de las columnas (por ejemplo, información sobre las muestras)
col_metadata <- data.frame(SampleID = paste0("Sample", 1:10), Condition = rep(c("Control", "Treatment"), 5))

# Crear el objeto SummarizedExperiment
se_object <- SummarizedExperiment(
  assays = list(counts = data_matrix),  # Asumimos que los datos son conteos
  rowData = row_metadata,               # Metadatos de las filas
  colData = col_metadata                # Metadatos de las columnas
)

# Ver el objeto creado
se_object





