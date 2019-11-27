

###########OBTENER INFORMACION DEL TCGA USANDO TCGABIOLINKS V 2.1.7
##########USAN R VANILLA 3.3.1
###CANCER DE CERVIX
##########19/ABRIL/2017

library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
query <- GDCquery(project = "TCGA-BRCA",
                  legacy = FALSE,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")

GDCdownload(query, method = "api")
data <- GDCprepare(query)

#Obtener todos los datos clinicos
#TCGA-BRCA

dataClin <- GDCquery_clinic(project = "TCGA-BRCA", "Clinical")

#####Extraer solo los read count

BRCAMatrix <- assay(data,"HTSeq - Counts")


#####VIEJO
#datos_read_count <- data[ , grepl( "read_count" , names( data ) ) ]
#Extraer los nombres de los mirnas de data y asignarlos a la extraccion de readcount
#rownames(datos_read_count) <- data$miRNA_ID

#Remover texto "readcounts del header de las columnas"

#nombres_mal <- colnames(datos_read_count)
#nombres_bien <- gsub("read_count_", "", nombres_mal)

#Reinsertar los nombres bien

#colnames(datos_read_count) <- nombres_bien

#write.csv(datos_read_count,"datos_read_count_BRCA.csv")

#----------------FIN DE OBTENCION DE LOS DATOS-----

#---------------Subset de los datos--------------------

#----------------Filtrado de tumores por tipo
#Que datos son PRIMARY SOLID TUMOR
#USAR ?TCGAquery_SampleTypes para ver codigos disponibles
datos_TP <- TCGAquery_SampleTypes(query$results[[1]]$cases,"TP") 
#Obtener la matriz
matriz_TP <- subset(BRCAMatrix, select=datos_TP)

write.csv(matriz_TP,"matriz_TP_BRCA.csv")

#Que datos son Solid Tissue Normal
datos_NT <- TCGAquery_SampleTypes(query$results[[1]]$cases,"NT") 
#Obtener la matriz
matriz_NT <- subset(BRCAMatrix, select=datos_NT)

write.csv(matriz_NT,"matriz_NT_BRCA.csv")

#--------------Filtrado mediante datos clinicos##

#########Matriz TP#####
#Obtencion de todos los barcodes de los TP
barcode_TP <- colnames(matriz_TP)

#Remover los ultimos 16 caracteres del barcode para que se ajusten a los datos clinicos
barcode_TP= substr(barcode_TP,1,nchar(barcode_TP)-16)

#Reemplazar los barcodes ajustados en la matriz

colnames(matriz_TP) <- barcode_TP
###############################

#Matriz NT############
#Obtencion de todos los barcodes de los NT
barcode_NT <- colnames(matriz_NT)

#Remover los ultimos 16 caracteres del barcode para que se ajusten a los datos clinicos
barcode_NT= substr(barcode_NT,1,nchar(barcode_NT)-16)

#Reemplazar los barcodes ajustados en la matriz

colnames(matriz_NT) <- barcode_NT
#####################################


#-----Organizacion de datos clinicos-------------------

#Obtener el bcr_patient_barcode,tumor_stage, vital_status,days_to_death,days_to_last_followup,progression_or_recurrence

#Remover columnas duplicadas de los datos clinicos
datos_clinicos_sinduplicados <- dataClin[ , !duplicated(colnames(dataClin))]

datos_clinicos_filtrados <- select(datos_clinicos_sinduplicados, bcr_patient_barcode,tumor_stage, vital_status, days_to_death,days_to_last_follow_up, progression_or_recurrence)

write.csv(datos_clinicos_filtrados,"datos_clinicos_filtrados_BRCA.csv")


##Extraccion de los nombres de los genes del summarized experiment
write.csv(matriz_NT,"matriz_NT_ENSEMBLE_BRCA_counts.csv")
write.csv(matriz_TP,"matriz_TP_ENSEMBLE_BRCA_counts.csv")

NOMBRES_GENES_MATRIZ <- rowData(data)$external_gene_name
rownames(matriz_NT) <- NOMBRES_GENES_MATRIZ
rownames(matriz_TP) <- NOMBRES_GENES_MATRIZ

#Guardar ambas matrices
write.csv(matriz_NT,"matriz_NT_GENES_BRCA_counts.csv")
write.csv(matriz_TP,"matriz_TP_GENES_BRCA_counts.csv")

