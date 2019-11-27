##########SCRIPT DATOS TCGA####
####    MODULO 3######
### Antonio Martinez
#Analisis de expresion diferencial
#Este script lee las dos matrices del modulo 2 y hace el analisis
#de expresion diferencial usando edgeR y DESeq2
rm(list=ls(all=TRUE))
#library(TCGAbiolinks)
library(DESeq2)
library("gplots")
library("Biobase")
library("RColorBrewer")

#Leer matrices generadas del modulo 2
matriz_filtro1 <- read.csv("matriz_NT_MATURE_MIRNAS_BRCA_readcount.csv",header = TRUE,row.names = 1)
matriz_filtro2 <- read.csv("matriz_stage_1_BRCA_readcount.csv",header = TRUE,row.names = 1)
matriz_filtro3 <- read.csv("matriz_stage_2_BRCA_readcount.csv",header = TRUE,row.names = 1)
matriz_filtro4 <- read.csv("matriz_stage_3_BRCA_readcount.csv",header = TRUE,row.names = 1)
matriz_filtro5 <- read.csv("matriz_stage_4_BRCA_readcount.csv",header = TRUE,row.names = 1)


nombre_filtro1 <- c("normales")
nombre_filtro2 <- c("tumores1")
nombre_filtro3 <- c("tumores2")
nombre_filtro4 <- c("tumores3")
nombre_filtro5 <- c("tumores4")



############Analisis con DESeq2#########################

#Contar el numero de pacientes en matriz filtrada uno y matriz filtrada 2

numero_pacientes_matriz1 <- ncol(matriz_filtro1)
numero_pacientes_matriz2 <- ncol(matriz_filtro2)
numero_pacientes_matriz3 <- ncol(matriz_filtro3)
numero_pacientes_matriz4 <- ncol(matriz_filtro4)
numero_pacientes_matriz5 <- ncol(matriz_filtro5)

#Unir matriz filtrada 1 con matriz filtrada 2

matriz_filtrada_unida <- cbind(matriz_filtro1,matriz_filtro2,
                               matriz_filtro3,matriz_filtro4,
                               matriz_filtro5)

#Guardar matriz filtrada unida

#write.csv(matriz_filtrada_unida,"matriz_filtrada_unida.csv")

#Leer matriz filtrada unida desde la variable

countdata <- matriz_filtrada_unida


# Convertir datos en matriz
countdata <- as.matrix(countdata)
#head(countdata)

# Asignar Condicion, diseño del experimento (primeros datos son el control, segundos el tratamiento)


(condition <- factor(c( rep(c(nombre_filtro1,nombre_filtro2,
                              nombre_filtro3,
                              nombre_filtro4,
                              nombre_filtro5), 
                            c(numero_pacientes_matriz1,
                              numero_pacientes_matriz2,
                              numero_pacientes_matriz3,
                              numero_pacientes_matriz4,
                              numero_pacientes_matriz5)))))

# Analisis con DESeq2 ----------------------------------------------------

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds


###vst transformation
vst <- varianceStabilizingTransformation(dds, blind=FALSE)
vst_valores <- assay(vst)
vst_grupo <- as.character(condition)

vst_union_grupo <- as.data.frame(rbind(vst_valores,vst_grupo))
write.csv(vst_union_grupo,"vst_transformation.csv")








