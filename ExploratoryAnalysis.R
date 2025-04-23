#You have your datasets feno and genes imported and tidied
#Now we focus on the analisis of this data

############################################
#Analysis of feno's stats
# Calculate the mean age and round the result
mean_age <- mean(feno$age, na.rm = TRUE)
round(mean_age)

# Display a summary of the age data (including min, median, mean, max, etc.)
summary(feno$age)

# Calculate the median age
median_age <- median(feno$age, na.rm = TRUE)
median_age

# Calculate the mean postoperative PSA and round the result
mean_psa_postop <- mean(feno$postop.psa, na.rm = TRUE)
round(mean_psa_postop)

# Calculate the mean preoperative PSA and round the result
mean_psa_preop <- mean(feno$preop.psa, na.rm = TRUE)
round(mean_psa_preop)



###########################################################################
# Specific questions to understand some characteristics of our patients


#Quantity of primary tumor tissue and normal prostatic tissue?
tipo_tejido <- table(feno$source_name_ch1)
tipo_tejido
# There is more primary tumor tissue than normal prostatic tissue

# How many patients of each age have each type of primary/secondary Gleason score?
table(feno$gleason.primary, feno$age)
table(feno$gleason.secondary, feno$age)
table(feno$gleason.sum, feno$age)
table(feno$biochemical.recurrence)
table(feno$tmprss2.erg)

# Age histogram
library(ggplot2)
feno$age

ggplot(feno, aes(x = age)) +
  geom_histogram(fill = "blue", color = "black") +
  labs(title = "Age Histogram", x = "Age", y = "Frequency") +
  theme_minimal()

ggsave("HistEdades.png", bg = "white", height = 4, width = 6, units = "in")

# Is there an association between the primary Gleason pattern and preoperative PSA? (jitter plot)
ggplot(feno, aes(x = gleason.primary, y = preop.psa)) +
  geom_jitter(width = 0.1)

# Is there an association between the secondary Gleason pattern and postoperative PSA? (jitter plot)
ggplot(feno, aes(x = gleason.secondary, y = postop.psa)) +
  geom_jitter(width = 0.1)

# Is there an association between the secondary Gleason pattern and preoperative PSA? (jitter plot)
ggplot(feno, aes(x = gleason.secondary, y = preop.psa)) +
  geom_jitter(width = 0.1)

# Is there an association between the primary Gleason pattern and postoperative PSA? (jitter plot)
ggplot(feno, aes(x = gleason.primary, y = postop.psa)) +
  geom_jitter(width = 0.1)




##########################################################################
library(skimr)
library(gtsummary)
library(tidyverse)
library(gt)
skim(feno)

# For getting a better summary of patient characteristics:
# General summary provided by R using table feno
tab1 <- feno %>% select(-title, -geo_accession, -gleason.primary, -gleason.secondary) %>% 
  tbl_summary()

gtsave(data = as_gt(tab1), filename = "prueba.gt.feno.png")
tab1




# Since I found irrelevant columns in my gt table, I decided to keep the relevant ones in a new list called lista_variables:

# See all column names in feno
colnames(feno)

# Create a list of all variables that contain relevant information from feno
variables <- colnames(feno)[c(3:9, 12:22)]
print(variables)

# Store these variables in a list
lista_variables <- as.list(variables)
print(lista_variables)




# Save phenotypic and genotypic tables as RDS files for later use

# Replace the path below with the folder where you want to save your data
saveRDS(feno, file = "path/to/your/folder/PhenotypicTable.rds")
saveRDS(genes, file = "path/to/your/folder/GenotypicTable.rds")






# Create a new table called feno1 from feno, keeping only the rows corresponding to RNA-Seq and primary tumor samples. So, pay attention in columns "library_strategy" and "source_name_ch1"). From now on, use this filtered data using these new tables.

feno1 <- feno %>%
  dplyr::filter(library_strategy == "RNA-Seq" & source_name_ch1 == "primary tumor")

# Adapt the genes table to match feno1; the result is called genes1
genes1 <- genes %>%
  dplyr::select('GeneID', any_of(feno1$geo_accession))



#Since multiple samples were collected from each patient, it was useful to create a new column called $ID.paciente, which extracts the patient ID from the first three digits of the feno$title column. This allowed us to confirm that the number of unique patients in the dataset matches the number reported in the original publication.


# Create a new column with patient ID, extracting the first 1 to 3 digits (optionally preceded by a letter)
feno$ID.paciente <- str_extract(feno$title, "^[A-Za-z]?[0-9]{1,3}")

# Create a second column with the remaining part of the title, representing the sample or analysis type
feno$ID.paciente.resto <- str_remove(feno$title, "^[A-Za-z]?[0-9]{1,3}")

# Reorder columns for easier visualization and sort by patient ID and sample type
feno <- feno %>%
  select(title, ID.paciente, ID.paciente.resto, source_name_ch1, library_strategy, everything()) %>%
  arrange(ID.paciente, source_name_ch1)



# Create a summary table of feno1 using gtsummary, excluding irrelevant columns
tab2 <- feno1 %>%
  select(-title, -geo_accession, -gleason.primary, -gleason.secondary, -ID.paciente, -ID.paciente.resto,
         -source_name_ch1, -library_strategy, -description, -library_selection) %>%
  tbl_summary()

# Save the summary table as a .png file
gtsave(data = as_gt(tab2), filename = "prueba.gt.feno1.png")


# Same gt summary table as before, but now stratified by a specific variable: biochemical.recurrence
tab3 <- feno1 %>%
  select(-title, -geo_accession, -gleason.primary, -gleason.secondary, -ID.paciente, -ID.paciente.resto,
         -source_name_ch1, -library_strategy, -description, -library_selection) %>%
  tbl_summary(by = biochemical.recurrence)

# Save the stratified summary table as a .png file
gtsave(data = as_gt(tab3), filename = "prueba.gt.biochrecu.png")
tab3


# Same gt summary table as before, but now stratified by a specific variable: risk.group
tab4=feno1 %>% select(-title,-geo_accession,-gleason.primary,-gleason.secondary, -ID.paciente, -ID.paciente.resto,-source_name_ch1, -library_strategy, -description, -library_selection) %>% 
  tbl_summary(by=risk.group)

# Save the stratified summary table as a .png file
gtsave(data = as_gt(tab4), filename = "prueba.gt.riskgroup.png")
tab4



saveRDS(feno1, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/TablaFenotipica1.rds")
saveRDS(genes1, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/TablaGenotipica1.rds")





# Boxplots to determine the age distribution of patients who experience recurrence

# 1. Evaluated using biochemical recurrence status
library(ggpubr)
ggplot(data = feno1 %>% filter(biochemical.recurrence != "NA"), 
       aes(x = biochemical.recurrence, y = age)) +
  geom_boxplot(fill = "violet") +
  labs(title = "Age distribution by biochemical recurrence status",
       x = "Biochemical recurrence", 
       y = "Age (years)") +
  theme_minimal() +
  scale_y_continuous(breaks = seq(floor(min(feno1$age, na.rm = TRUE) / 5) * 5, 
                                  ceiling(max(feno1$age, na.rm = TRUE) / 5) * 5, 
                                  by = 5)) +
  stat_compare_means(method = "t.test", label = "p.format")
ggsave("boxplot_age_biochrecu.png", width = 8, height = 6, dpi = 300)

# 2. Evaluated based on postoperative PSA levels
# Need to determine a PSA threshold that is clinically relevant for recurrence.
# For now, we use an arbitrary cutoff of 0.30
ggplot(feno1, aes(x = postop.psa >= 0.3, y = age)) +
  geom_boxplot(fill = "violet") +
  scale_x_discrete(labels = c(">0.3", "<=0.3")) +
  labs(title = "Age distribution by postoperative PSA levels",
       x = "Postoperative PSA levels", 
       y = "Age (years)") +
  theme_minimal()
ggsave("boxplot_edades_postopPSA.png", width = 8, height = 6, dpi = 300)

# 3. Evaluated based on residual tumor status (R0 vs R1)
ggplot(feno1, aes(x = r.residual_tumor, y = age)) +
  geom_boxplot(fill = "violet") +
  labs(title = "Age distribution by residual tumor status", 
       x = "Residual tumor", 
       y = "Age (years)") +
  theme_minimal()
ggsave("boxplot_edades_residualtumor.png", width = 8, height = 6, dpi = 300)

# Is there an association between age and preoperative PSA? (scatter plot)
library(ggplot2)

ggplot(feno1, aes(x = age, y = preop.psa, color = gleason.sum)) +
  geom_jitter(alpha = 0.6, width = 0.2, height = 0.2) +
  geom_smooth(method = "lm", color = "black", se = TRUE, fill = "gray70") + 
  labs(title = "Relationship between age and preoperative PSA levels",
       x = "Age (years)", 
       y = "Preoperative PSA levels",
       color = "Gleason Score") +
  theme_minimal() +
  scale_color_brewer(palette = "Paired") 

ggsave("jitterplot_edades_PSApreoperatorio.png", width = 8, height = 6, dpi = 300)





#Logistic regresion

# This step is FUNDAMENTAL because it defines "brf" as improvement and "bcr" as deterioration
table(feno1$biochemical.recurrence)
class(feno1$biochemical.recurrence)
feno1$biochemical.recurrence <- factor(feno1$biochemical.recurrence, levels = c("brf", "bcr"))


# Logistic regression between biochemical recurrence and preoperative PSA
modelo <- glm(biochemical.recurrence ~ preop.psa, data = feno1, family = binomial)
summary(modelo)

# Logistic regression between biochemical recurrence and Gleason score
modelo2 <- glm(biochemical.recurrence ~ gleason.sum, data = feno1, family = binomial)
summary(modelo2)



#To perform a logistic regression between gene expression levels and the variable biochemical.recurrence, we first need to combine the phenotypic and genotypic datasets.
#We start by selecting only the genes of interest: AR, YWHAZ, NDRG1, APOE, CRIP2, and


# Selected genes by their GeneID
muestra <- c(367, 7534, 10397, 348, 1397, 10631)

# Find indices of the genes of interest
indices <- which(genes1$GeneID %in% muestra)

# Subset the gene expression table using those indices
nueva_tabla <- genes1[indices,]
View(nueva_tabla)

# Rename gene IDs to gene names
nueva_tabla[1,1] <- "YWHAZ"
nueva_tabla[2,1] <- "NDRG1"
nueva_tabla[3,1] <- "POSTN"
nueva_tabla[4,1] <- "CRIP2"
nueva_tabla[5,1] <- "APOE"
nueva_tabla[6,1] <- "AR"

# Set gene names as row names and remove GeneID column
row.names(nueva_tabla) <- nueva_tabla$GeneID
nueva_tabla = nueva_tabla[-1]

# Transpose to have patients in rows and genes in columns
genes_translocados <- t(nueva_tabla)
View(genes_translocados)

# Save the transposed table
write.csv(genes_translocados, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/genes_translocados.csv", row.names = TRUE)
saveRDS(genes_translocados, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/nueva_tabla_translocada.rds")



#merge the clinical and molecular data using `merge`. I use `geo_accession` as the reference column that contains the patient IDs.
tabla_unificada <- merge.data.frame(feno1, genes_translocados, by.x = "geo_accession", by.y = "row.names", all = TRUE)
View(tabla_unificada)

#save the unified table as an RDS file
saveRDS(tabla_unificada, "C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/tabla_unificada.rds")




###############################################################################
#Evaluo regresion logistica entre genes y biochemical.recurrence
#linea 311 de analisisdatocrudosmolecularesyclinicosGSE229904
