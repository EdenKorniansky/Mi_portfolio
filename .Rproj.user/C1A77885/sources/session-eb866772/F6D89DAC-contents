setwd("C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904")

# import datos expresion genetica
genes <- read.delim("C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/GSE229904_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
View(genes)

# import datos fenotipicos
feno <- read.csv("C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/GSE229904_Annotations.csv", row.names=1)
View(feno)

#elimino columnas que no necesito
library(tidyverse)
feno <- feno %>%
  dplyr::select(!starts_with("character"))

#gsub("\\.ch1$", "", colnames(tabla)): Reemplaza el sufijo .ch1 al final ($ indica "al final de la cadena") con una cadena vacía ("").
colnames(feno) <- gsub("\\.ch1$", "", colnames(feno))

#factores de cada columna de gleason
cols_to_factor <- c("gleason.primary", "gleason.secondary", "gleason.sum", "risk.group",  
                    "source_name_ch1", "description", "library_selection", "library_strategy",  
                    "biochemical.recurrence", "cm", "isup", "lapc", "pn", "pt",  
                    "r.residual_tumor", "tmprss2.erg")

feno[cols_to_factor] <- lapply(feno[cols_to_factor], factor)

#corroboro si esta en formato numerico
feno$age
feno$gleason.primary
feno$gleason.secondary
feno$gleason.sum
feno$isup
feno$ln.percent
feno$postop.psa
feno$preop.psa


#corroboro si hay columnas que tienen la misma informacion y por lo tanto, deberia borrarlas
table(feno$description)
table(feno$library_selection)
table(feno$library_strategy)
table(feno$cm)
table(feno$lapc)
table(feno$r.residual_tumor)
table(feno$risk.group)


#promedios y medianas de edad y psa
prom_edad <- mean(feno$age, na.rm = T)
round(prom_edad)
#forma mas eficiente de verlo
summary(feno$age)

mediana_edad <- median(feno$age, na.rm = T)
mediana_edad

prom_psa_postop <- mean(feno$postop.psa, na.rm = T)
round(prom_psa_postop)

prom_psa_preop <- mean(feno$preop.psa, na.rm = T)
round(prom_psa_preop)



#Hay mayor tejido de tumor primario que de tejido normal prostatico
tipo_tejido <- table(feno$source_name_ch1)
tipo_tejido

#cuantos pacientes de cada edad tienen gleason primario/secundario de cada tipo
table(feno$gleason.primary, feno$age)
table(feno$gleason.secondary, feno$age)
table(feno$gleason.sum, feno$age)
table(feno$biochemical.recurrence)
table(feno$tmprss2.erg)
#Histograma de edad
library(ggplot2)
feno$age

ggplot(feno, aes(x = age)) +
  geom_histogram(fill = "blue", color = "black") +
  labs(title = "Histograma de edades", x = "edades", y = "Frecuencia") +
  theme_minimal()

ggsave("HistEdades.png", bg = "white", height = 4, width = 6, units = "in")

#hay una asosiacion entre el primary psa pattern y el PSA preoperatorio? (jitter plot)

ggplot(feno, aes(x=gleason.primary, y=preop.psa))+
  geom_jitter(width = 0.1)

#hay una asosiacion entre el secondary psa pattern y el PSA postperatorio? (jitter plot)

ggplot(feno, aes(x=gleason.secondary, y=postop.psa))+
  geom_jitter(width = 0.1)

#hay una asosiacion entre el secondary psa pattern y el PSA preoperatorio? (jitter plot)

ggplot(feno, aes(x=gleason.secondary, y=preop.psa))+
  geom_jitter(width = 0.1)




#hay una asosiacion entre el primary psa pattern y el PSA postperatorio? (jitter plot)

ggplot(feno, aes(x=gleason.primary, y=postop.psa))+
  geom_jitter(width = 0.1)
library(skimr)
library(gtsummary)
library(tidyverse)
library(gt)
skim(feno)


#Resumen general que me brinda R
tab1=feno %>% select(-title,-geo_accession,-gleason.primary,-gleason.secondary) %>% 
  tbl_summary()

gtsave(data = as_gt(tab1), filename = "prueba.gt.feno.png")
tab1


# para ver las variables que tengo
colnames(feno)


#una lista de todas las variables que tengan informacion relevante de info feno
variables <- colnames(feno)[c(3:9, 12:22)]
print(variables)

#luego las inserto en una lista
lista_variables <- as.list(variables)
print(lista_variables)


saveRDS(feno, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/TablaFenotipica.rds")
saveRDS(genes, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/TablaGenotipica.rds")



#creo columna de id paciente queandome con los primeros 3 digitos
feno$ID.paciente <- str_extract(feno$title, "^[A-Za-z]?[0-9]{1,3}")

#colunda del resto, hace lo mismo pero removiendo los 3 digitos para quedarse con el tipo de analisis de la muestra
feno$ID.paciente.resto<- str_remove(feno$title, "^[A-Za-z]?[0-9]{1,3}")

# solo para ver todo juntos, ordena las columnas y ordena en modo ascendente el ID paciente junto con el tejido extraido
feno=feno %>% 
  select(title,ID.paciente, ID.paciente.resto, source_name_ch1, library_strategy, everything()) %>% 
  arrange(ID.paciente, source_name_ch1)



#Creo una tabla llamada feno1 a partir de feno, esta tiene que dejar aquellas filas que sean RNAseq y tumor primario

feno1 <- feno %>%
  dplyr::filter(library_strategy == "RNA-Seq" & source_name_ch1 == "primary tumor")

#la tabla genes se amolda a feno1, la llamamos genes1
genes1 <- genes %>%
  dplyr::select('GeneID', any_of(feno1$geo_accession))

#Resumen con tabla gt de feno1
tab2=feno1 %>% select(-title,-geo_accession,-gleason.primary,-gleason.secondary, -ID.paciente, -ID.paciente.resto,-source_name_ch1, -library_strategy, -description, -library_selection) %>% tbl_summary()

gtsave(data = as_gt(tab2), filename = "prueba.gt.feno1.png")



#tabla summary estratificada por biochemical recurrance divido en dos poblaciones: bcr y brf
tab3=feno1 %>% select(-title,-geo_accession,-gleason.primary,-gleason.secondary, -ID.paciente, -ID.paciente.resto,-source_name_ch1, -library_strategy, -description, -library_selection) %>% 
  tbl_summary(by=biochemical.recurrence)


gtsave(data = as_gt(tab3), filename = "prueba.gt.biochrecu.png")
tab3

# tabla summary estratificada por risk group
tab4=feno1 %>% select(-title,-geo_accession,-gleason.primary,-gleason.secondary, -ID.paciente, -ID.paciente.resto,-source_name_ch1, -library_strategy, -description, -library_selection) %>% 
  tbl_summary(by=risk.group)


gtsave(data = as_gt(tab4), filename = "prueba.gt.riskgroup.png")
tab4


saveRDS(feno1, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/TablaFenotipica1.rds")
saveRDS(genes1, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/TablaGenotipica1.rds")
feno1 <- feno1 %>% filter(!is.na(age))



#boxplots para determinar cual es la edad de aquellos pacientes que recaen
#1 evaluado con biochemical.recurrance
library(ggpubr)
ggplot(data = feno1 %>% filter(biochemical.recurrence != "NA"), 
       aes(x = biochemical.recurrence, y = age)) +
  geom_boxplot(fill = "violet") +
  labs(title = "Edades en función de BCR y BFR",
       x = "Biochemical recurrence", 
       y = "Age (years)") +
  theme_minimal() +
  scale_y_continuous(breaks = seq(floor(min(feno1$age, na.rm = TRUE) / 5) * 5, 
                                  ceiling(max(feno1$age, na.rm = TRUE) / 5) * 5, 
                                  by = 5)) +
  stat_compare_means(method = "t.test", label = "p.format")
ggsave("boxplot_age_biochrecu.png", width = 8, height = 6, dpi = 300)

#2 evaluado en psa postoperatorio 
# Averiguar que valor de PSA postoperatorio seria alto y tenga relacion con recaida
#Por ahora uso un ejemplo de 0,30
ggplot(feno1, aes(x = postop.psa >= 0.3, y = age)) +
  geom_boxplot(fill = "violet") +
  scale_x_discrete(labels = c(">0.3", "<=0.3")) +
  labs(title = "Edades en función de niveles de PSA postoperatorio",
       x = "Niveles de PSA postoperatorio", 
       y = "Edades") +
  theme_minimal()
ggsave("boxplot_edades_postopPSA.png", width = 8, height = 6, dpi = 300)

#3 Evaluado con r.residual_tumor
ggplot(feno1, aes(x = r.residual_tumor, y = age)) +
  geom_boxplot(fill = "violet") +
  labs(title = "Edades en funcion de R0 y R1",x = "Residual tumor", y = "Edades") +
  theme_minimal()

ggsave("boxplot_edades_residualtumor.png", width = 8, height = 6, dpi = 300)



#Hay una asosiacion entre la edad y el PSA preoperatorio? (esto se hace con un grafico de puntos)
library(ggplot2)

ggplot(feno1, aes(x = age, y = preop.psa, color = gleason.sum)) +
  geom_jitter(alpha = 0.6, width = 0.2, height = 0.2) +
  geom_smooth(method = "lm", color = "black", se = TRUE, fill = "gray70") + 
  labs(title = "Edad en función de PSA preoperatorio",
       x = "Age (years)", 
       y = "Preoperative PSA levels",
       color = "Gleason Score") +  # Etiqueta de la leyenda
  theme_minimal() +
  scale_color_brewer(palette = "Paired") 


ggsave("jitterplot_edades_PSApreoperatorio.png", width = 8, height = 6, dpi = 300)




#FUNDAMENTAL ESTE PASO
table(feno1$biochemical.recurrence)
class(feno1$biochemical.recurrence)
feno1$biochemical.recurrence <- factor(feno1$biochemical.recurrence, levels = c("brf", "bcr"))

#Regresion logistica entre biochemical.recurrence y preop.psa
modelo <- glm(biochemical.recurrence~preop.psa, data = feno1, family = binomial)
summary(modelo)

#Regresion logistica entre biochemical.recurrence y gleason.sum
modelo2 <- glm(biochemical.recurrence~gleason.sum, data = feno1, family = binomial)
summary(modelo2)


#Creacion de tabla unificada de feno1 y genes(AR, YWHAZ, NDRG1, APOE, CRIP2, POSTN)

muestra <- c(367, 7534, 10397, 348, 1397, 10631)

indices <- which(genes1$GeneID %in% muestra)

#uso los indices para crear la tabla nueva
nueva_tabla<-genes1[indices,]
View(nueva_tabla)

#cambio los numeros de ID de los genes por sus nombres en la tabla nueva
nueva_tabla[1,1] <- "YWHAZ"
nueva_tabla[2,1] <- "NDRG1"
nueva_tabla[3,1] <- "POSTN"
nueva_tabla[4,1] <- "CRIP2"
nueva_tabla[5,1] <- "APOE"
nueva_tabla[6,1] <- "AR"

#cambio el nombre de las columnas de los extremos por los nombres de los genes
row.names(nueva_tabla) <- nueva_tabla$GeneID
nueva_tabla = nueva_tabla[-1]

genes_translocados <- t(nueva_tabla)
View(genes_translocados)
write.csv(genes_translocados, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/genes_translocados.csv", row.names = TRUE)
saveRDS(genes_translocados, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/nueva_tabla_translocada.rds")


#unifico la clinica con la molecular usando merge, uso geo_accession como referencia de donde arranca la columna que contiene a los pacientes
tabla_unificada <- merge.data.frame(feno1, genes_translocados, by.x = "geo_accession", by.y = "row.names", all = TRUE)
View(tabla_unificada)

#guardo la tabla unificada en un archivo RDS
saveRDS(tabla_unificada,"C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/tabla_unificada.rds")


# I merge the clinical and molecular data using `merge`. I use `geo_accession` as the reference column that contains the patient IDs.
tabla_unificada <- merge.data.frame(feno1, genes_translocados, by.x = "geo_accession", by.y = "row.names", all = TRUE)
View(tabla_unificada)

# I save the unified table as an RDS file
saveRDS(tabla_unificada, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/tabla_unificada.rds")





###############################################################################
#Evaluo regresion logistica entre genes y biochemical.recurrence

library(stats)
# Seleccionar los nombres de los genes
genes <- colnames(tabla_unificada)[25:30]
genes
# Inicializar una lista para guardar los modelos
modelos <- list()
print(colnames(tabla_unificada))

# Loop para ajustar modelos de regresión logística para cada gen individualmente
for (gen in genes) {
  formula <- as.formula(paste("biochemical.recurrence ~", gen))  # Convertir string en fórmula
  modelos[[gen]] <- glm(formula, data = tabla_unificada, family = binomial)
  
  # Mostrar resumen del modelo
  print(summary(modelos[[gen]]))
}




#realizar una tabla con 3 columnas: genes, estimate y p-value

# Inicializar una lista vacía para almacenar los resultados
resultados <- list()

# Iterar sobre los genes y extraer coeficiente y p-valor
for (gen in genes) {
  modelo <- modelos[[gen]]
  resumen <- summary(modelo)
  
  # Extraer el coeficiente (Estimate) y el p-valor (Pr>|z|)
  coeficiente <- coef(resumen)[gen, "Estimate"]
  p_valor <- coef(resumen)[gen, "Pr(>|z|)"]
  
  # Guardar en la lista
  resultados[[gen]] <- data.frame(Gene = gen, Estimate = coeficiente, P_Value = p_valor)
}

# Convertir la lista en un data frame final
tabla_gen_est_pv<- do.call(rbind, resultados) #Esta convirtiendo una lista de data frames (llamada resultados) en un solo data frame
rownames(tabla_gen_est_pv) <- NULL  # Elimina nombres de fila si están generando problemas


#Grafico
# Cargar la librería ggplot2
library(ggplot2)

# Crear el gráfico

ggplot(tabla_gen_est_pv, aes(x = Estimate, y = reorder(Gene, Estimate), fill = P_Value)) +
  geom_col(width = 0.4, color = "black") +  # Ajustar el ancho de las barras (0.5 es más delgado)
  scale_fill_gradient(low = "black", high = "white", name = "P-Value") +  # Gradiente de color
  scale_x_continuous(limits = c(-0.05, 0.05))+
  theme_minimal() +  # Tema limpio
  labs(x = "Estimate", y = NULL, title = "Coeficientes de Regresión y su Significancia") +
  theme(axis.text.y = element_blank()) +  # Eliminar etiquetas del eje Y
  geom_text(aes(label = Gene), hjust = ifelse(tabla_gen_est_pv$Estimate > 0, -0.2, 1.2), size = 4) +  # Nombres al lado de las barras
  theme(legend.position = "right")  # Ubicar la leyenda a la derecha


ggsave("tablaGenesEstimatePvalue.png",width = 7, height = 5, bg = "white")





###################################################################################
#Regresion logistica entre biochemical.recurrence y datos clinicos + graficos

#inicializo variable datos con los datos clinicos QUE ME INTERESAN
datos <- c("age", "cm", "gleason.sum", "isup", "lapc","ln.percent","pn", "postop.psa", "preop.psa", "pt", "r.residual_tumor", "risk.group","tmprss2.erg")
library(stats)


# Inicializar la lista de modelos
modelos_clinicos <- list()

# Loop para ajustar modelos de regresión logística para cada dato individualmente
for (dato in datos) {
  formula <- as.formula(paste("biochemical.recurrence ~", dato))  # Convertir string en fórmula
  modelos_clinicos[[dato]] <- glm(formula, data = tabla_unificada, family = binomial)
  
  # Mostrar resumen del modelo
  print(summary(modelos_clinicos[[dato]]))
}



# Inicializar una lista vacía para almacenar los resultados
# Crear una lista para almacenar los resultados
resultados <- list()

# Iterar sobre los datos y extraer coeficiente y p-valor
for (dato in datos) {
  modelo <- modelos_clinicos[[dato]]
  resumen <- summary(modelo)
  coeficientes <- coef(resumen)  # Extraer la tabla de coeficientes
  
  # Buscar el nombre correcto del coeficiente (manejo de nombres inconsistentes)
  coef_nombres <- rownames(coeficientes)
  coef_indice <- grep(paste0("^", dato), coef_nombres)  # Encuentra la fila correcta
  
  if (length(coef_indice) > 0) {  # Verifica si encontró el coeficiente
    coeficiente <- coeficientes[coef_indice[1], "Estimate"]
    p_valor <- coeficientes[coef_indice[1], "Pr(>|z|)"]
    
    # Guardar en la lista
    resultados[[dato]] <- data.frame(Dato = dato, Estimate = coeficiente, P_Value = p_valor)
  } else {
    cat("⚠️ Advertencia: No se encontró coeficiente para", dato, "\n")
  }
}

# Convertir la lista en un data frame final
tabla_final <- do.call(rbind, resultados)
rownames(tabla_final) <- NULL  # Elimina nombres de fila si están generando problemas
view(tabla_final)


# Gráfico
library(ggplot2)

ggplot(tabla_final, aes(x = Estimate, y = reorder(Dato, Estimate), fill = P_Value)) +
  geom_col(width = 0.6, color = "black") +  # Agregar borde negro para mayor visibilidad
  scale_x_continuous(expand = c(0, 0)) + # asegura que la escala siempre este ajustada a los valores en X
  scale_fill_gradient(low = "black", high = "white", name = "P-Value") +  # Gradiente de color
  theme_minimal() +  # Tema limpio
  labs(x = "Estimate", y = NULL, title = "Coeficientes de Regresión de Datos Clínicos") +
  theme(axis.text.y = element_text(size = 10)) +  # Restaurar etiquetas del eje Y 
  theme(legend.position = "right")  # Ubicar la leyenda a la derecha


# Guardar el gráfico
ggsave(filename = "RegresionLogistica_BRC_DatosClinicos.png", width = 7, height = 5, bg = "white")




#Grafico igual al anterior pero que se note bien la diferencia de p-value con colores

#Grafico con p-value >= 0,05 blanco y <0,05 color oscuro/claro
#p >= 0,05, Lo asignamos como NA, y luego usamos na.value = "white" en la escala de color.
#p < 0,05, Usamos el propio valor de P_Value para el color
#xq scale_fill_gradient(), los valores deben ser numéricos

tabla_final$color <- ifelse(tabla_final$P_Value >= 0.05, NA, tabla_final$P_Value)

ggplot(tabla_final, aes(x = Estimate, y = reorder(Dato, Estimate), fill = color)) +
  geom_col(width = 0.6, color = "black") +  # Borde negro para todas las barras
  scale_x_continuous(expand = c(0, 0)) +  # Ajustar escala en X
  scale_fill_gradient(low = "pink", high = "red", name = "P-Value", na.value = "white") +  
  theme_minimal() +  
  labs(x = "Estimate", y = NULL, title = "Coeficientes de Regresión de Datos Clínicos") +
  theme(axis.text.y = element_text(size = 10), legend.position = "right")


# Guardar el gráfico
ggsave(filename = "RegresionLogistica_BRC_DatosClinicosConColores.png", width = 14, height = 5, bg = "white")




##############################################################################
#Analisis de expresion diferencial --> differential gene expression DGE
#Analiso biochemical recurrence frente a toda la expresion de mis genes



#Lo que debo hacer es utilizar la tabla genes1 que queda estatica  y feno1 translocada con solo la variable geo_accession


#tengo en las filas los nombres de genes y columnas nombres de muestras
genesDGE <- genes1
View(genesDGE)
rownames(genesDGE) <- as.character(genesDGE$GeneID)  #nombres de los genes
genesDGE <- genesDGE[,-1]


#Me quedo solo con los nombres de las muestras y biochemical recurrence
fenoDGE <- data.frame(muestra = rownames(feno1), biochemical.recurrence = feno1$biochemical.recurrence)
View(fenoDGE)

saveRDS(fenoDGE, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/fenoDGE.rds")
saveRDS(genesDGE, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/genesDGE.rds")

# necesitamos una matriz de expresion que se va a llamar count_norm (ya normalizada) y un vector que sea la condicion/grupo al que pertenece cada paciente, que se va a llamar conditions, que sea de tipo factor

count_norm <- genesDGE
conditions <- fenoDGE$biochemical.recurrence
  
# asegurarnos de que count_norm sea un dataframe
count_norm <- as.data.frame(count_norm)


# Run the Wilcoxon rank-sum test for each gene
pvalues <- sapply(1:nrow(count_norm), function(i){
  data <- cbind.data.frame(gene = as.numeric(t(count_norm[i,])), conditions)
  p <- wilcox.test(gene~conditions, data)$p.value
  return(p)
})
#cbind.data.frame(...): crea un data frame con dos columnas:
#gene: expresión del gen
#conditions: condición de cada muestra



#ajuste común cuando se realizan múltiples pruebas estadísticas para reducir el riesgo de obtener falsos positivos
fdr <- p.adjust(pvalues, method = "fdr")


# Calculate the fold-change for each gene

conditionsLevel <- levels(conditions)
dataConBRF <- count_norm[,c(which(conditions==conditionsLevel[1]))]
dataConBCR <- count_norm[,c(which(conditions==conditionsLevel[2]))]

# se calcula el fold-change como el cociente de las medias
foldChanges <- log2(rowMeans(dataConBCR)/rowMeans(dataConBRF))
#Calcula el log2 del fold-change gen a gen
#Si el resultado es >0, el gen está sobreexpresado en BCR.
#Si es <0, el gen está subexpresado en BCR.

# Output results based on the FDR threshold 0.05
outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)
rownames(outRst) <- rownames(count_norm)
outRst <- na.omit(outRst)
view(outRst)

saveRDS(outRst, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/outRst.rds")


#cuantos genes hay con un FDR menor a 0.05?

# Genes con FDR < 0.05
significativos_005 <- which(fdr < 0.05, arr.ind = TRUE)
cantidad_significativos_005 <- length(significativos_005) 
cantidad_significativos_005 #1253


#cuales son los 20 genes con el menor FDR? --> los mas significativos!!!
orden_fdr <- order(fdr, na.last = NA)

veinte_menores_FDR <- orden_fdr[1:20]
veinte_menores_FDR



#cuales son los 20 genes con el mayor valor absoluto de logfoldchange?
orden_logFC <- order(abs(foldChanges), decreasing = TRUE, na.last = NA)
veinte_mayores_LFC <- orden_logFC[1:20]
veinte_mayores_LFC


#resumir estos datos: volcanoplot

library(ggplot2)
library(ggrepel)  # Para agregar etiquetas sin superposición

# Crear una columna para clasificar genes según FDR y logFC
outRst$Significativo <- ifelse(outRst$FDR < 0.05 & abs(outRst$log2foldChange) > 1, 
                               "Significativo", "No significativo")


# Obtener los 20 genes más significativos, osea, con menor FDR
top_genes <- rownames(outRst)[orden_fdr[1:20]]

# Volcano plot
outRst$gene <- rownames(outRst)

volcanoplot <- ggplot(outRst, aes(x = log2foldChange,
                                  y = -log10(FDR),
                                  color = Significativo,
                                  text = paste0("Gen: ", gene,
                                                "\nlog2FC: ", round(log2foldChange, 2),
                                                "\nFDR: ", signif(FDR, 3)))) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Significativo" = "red", "No significativo" = "gray")) +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-log10(FDR)") +
  theme_minimal() +
  geom_text_repel(data = outRst[top_genes, ],
                  aes(label = gene),
                  size = 3, max.overlaps = 10)

ggsave("C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/volcano_plot.png", plot = volcanoplot, width = 8, height = 6, dpi = 300)


library(plotly)
ggplotly(volcanoplot, tooltip = "text")
