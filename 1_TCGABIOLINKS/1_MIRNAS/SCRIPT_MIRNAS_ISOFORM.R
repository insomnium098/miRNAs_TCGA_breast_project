

###########OBTENER INFORMACION DEL TCGA USANDO TCGABIOLINKS
##########USAN R VANILLA 3.3.1
###BRCA
library(SummarizedExperiment)
library(dplyr)
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-BRCA",
                  legacy = TRUE,
                  data.category = "Gene expression",
                  data.type = "miRNA isoform quantification")


GDCdownload(query, method = "client")
data <- GDCprepare(query)

#Obtener todos los datos clinicos
#TCGA-BRCA

dataClin <- GDCquery_clinic(project = "TCGA-BRCA", "Clinical")

#####Extraer solo los read count

matriz_readcount_t <- data #assay(data,"HTSeq - FPKM")


##########################3SCRIPT PROCESAMIENTO DATOS MIRNAS
#ESTE SCRIPT PROCESA EL DATAFRAMDE DEL TCGABIOLINKS Y PERMITE CONTINUAR CON EL PIPELINE

#TRANSPONER MATRIZ
matriz_readcount_t_TRANSPOSED <- t(matriz_readcount_t)

#extraer nombres de los mirnas

nombres_mirnas <- matriz_readcount_t$miRNA_ID

#EXTRAER LOS READ COUNT

index_matriz_readcount <- grep("read_count_",rownames(matriz_readcount_t_TRANSPOSED))

matriz_readcount <- matriz_readcount_t_TRANSPOSED[index_matriz_readcount,]

#transpose matriz_readcount

matriz_readcount_t <- t(matriz_readcount)

#agregar rownames

rownames(matriz_readcount_t) <- nombres_mirnas

#REMOVER TEXTO READS PER MILLION

nombres_mal <- colnames(matriz_readcount_t)

nombres_bien <- gsub("read_count_","",nombres_mal,fixed = TRUE)

colnames(matriz_readcount_t) <- nombres_bien

###################################FIN DE PROCESAMIENTO DE MIRNAS

 
#Obtener la matriz
matriz_TP <- matriz_readcount_t#subset(matriz_readcount_t, select=datos_TP)

write.csv(matriz_TP,"matriz_TP_BRCA_SIN_16.csv")

#Que datos son Solid Tissue Normal

datos_NT <- TCGAquery_SampleTypes(query$results[[1]]$cases,"NT") 
#Obtener la matriz
matriz_NT <- subset(matriz_readcount_t, select=datos_NT)

write.csv(matriz_NT,"matriz_NT_BRCA_MIRNAS_SIN_16.csv")


#-----Organizacion de datos clinicos-------------------

#Obtener el bcr_patient_barcode,tumor_stage, vital_status,days_to_death,days_to_last_followup,progression_or_recurrence

#Remover columnas duplicadas de los datos clinicos
datos_clinicos_sinduplicados <- dataClin[ , !duplicated(colnames(dataClin))]

datos_clinicos_filtrados <- select(datos_clinicos_sinduplicados, bcr_patient_barcode,tumor_stage, vital_status, days_to_death,days_to_last_follow_up, progression_or_recurrence)

write.csv(datos_clinicos_filtrados,"datos_clinicos_filtrados_BRCA.csv")


