# evaluando porque hay tantas muestras

TablaFenotipica <- readRDS("C:/Users/edenk/OneDrive/Pasantía Cáncer de Próstata/GSE229904/TablaFenotipica.rds")
datos=TablaFenotipica # solo para trabajar con una tabla que se llame datos, que es mas rapido

# paquetes necesarios
library(skimr)
library(tidyverse)


skim(datos)

# aprece ser que la columna title tiene info sobre el paciente del cual se toma la muestra. Evidentemente el nro que aparece al principio de cada celda, es el dentificador unico de cada paciente, y de cada paciente se tienen varias muetras, algunas de tumor, algnas de normal, y a su vez, a esas meustras se las analiza o por TNAseq o por miRNA seq


#[A-Za-z]? → Busca cero o una letra (opcional) al principio. Puede ser mayúscula o minúscula.
#[0-9]{1,3} → Busca entre 1 y 3 dígitos numéricos consecutivos.

datos$title <- gsub("^U", "", datos$title)
datos$ID.paciente <- str_extract(datos$title, "^[A-Za-z]?[0-9]{1,3}")



#hace lo mismo pero removiendo los 3 digitos para quedarse con el tipo de analisis de la muestra
datos$ID.paciente.resto<- str_remove(datos$title, "^[A-Za-z]?[0-9]{1,3}")

# solo para ver todo juntos, ordena las columnas y ordena en modo ascendente el ID paciente junto con el tejido extraido
datos=datos %>% 
  select(title,ID.paciente, ID.paciente.resto, source_name_ch1, library_strategy, everything()) %>% 
  arrange(ID.paciente, source_name_ch1)


# parece que hay 158 pacientes
unique(datos$ID.paciente) %>% length()


# hay 22 que empiezan con U, no se que significa eso
grepl(pattern = "u", ignore.case = T, x = datos$ID.paciente ) %>% sum()

# cuantas veces aparece cada una?
sort(table(datos$ID.paciente), decreasing = T)



cantidad_pacientes <- nrow(subset(datos, pn == "N0" & lapc == "lapc"))
cantidad_pacientes
