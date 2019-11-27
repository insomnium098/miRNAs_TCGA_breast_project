##########SCRIPT DATOS TCGA####
####    MODULO 3######
### Antonio Martinez
#Analisis de expresion diferencial
#Este script lee las dos matrices del modulo 2 y hace el analisis
#de expresion diferencial usando DESeq2
rm(list=ls(all=TRUE))
#library(TCGAbiolinks)
library(DESeq2)
library("gplots")
library("Biobase")
library("RColorBrewer")

#Leer matrices generadas del modulo 2
matriz_filtro1 <- read.csv("matriz_NT_MATURE_MIRNAS_BRCA_readcount.csv",header = TRUE,row.names = 1)
matriz_filtro2 <- read.csv("matriz_stage_1_BRCA_readcount.csv",header = TRUE,row.names = 1)

nombre_filtro1 <- c("normales")
nombre_filtro2 <- c("tumores")



############Analisis con DESeq2#########################

#Contar el numero de pacientes en matriz filtrada uno y matriz filtrada 2

numero_pacientes_matriz1 <- ncol(matriz_filtro1)
numero_pacientes_matriz2 <- ncol(matriz_filtro2)

#Unir matriz filtrada 1 con matriz filtrada 2

matriz_filtrada_unida <- cbind(matriz_filtro1,matriz_filtro2)

#Guardar matriz filtrada unida

#write.csv(matriz_filtrada_unida,"matriz_filtrada_unida.csv")

#Leer matriz filtrada unida desde la variable

countdata <- matriz_filtrada_unida


# Convertir datos en matriz
countdata <- as.matrix(countdata)
#head(countdata)

# Asignar Condicion, diseño del experimento (primeros datos son el control, segundos el tratamiento)


(condition <- factor(c( rep(c(nombre_filtro1,nombre_filtro2), c(numero_pacientes_matriz1,numero_pacientes_matriz2)))))

# Analisis con DESeq2 ----------------------------------------------------

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)#, parallel = TRUE)

# Plot dispersions
#png("qc-dispersions.png", 1000, 1000, pointsize=20)
#plotDispEsts(dds, main="Dispersion plot")
#dev.off()

# Regularizar la transformacion de log en rlog para poder graficar en heatmap
#rld <- rlogTransformation(dds)
# Extraer la transformacion del s4class rld
#rld3 <- as.matrix(t(assay(rld)))
#Guardar transformacion en rld.csv
#cambiar texto "." por "-" del header de las columnas de la matriz rld3

#nombres_mal_rld3 <- row.names(rld3)
#nombres_bien_rld3 <- gsub(".", "-", nombres_mal_rld3,fixed = TRUE)
#row.names(rld3) <- nombres_bien_rld3
#write.csv(rld3,"rld.csv")
#head(assay(rld))
#hist(assay(rld))

# Obtener los datos de Differential Expression
res <- results(dds)
table(res$padj<0.05)
## Ordenar por adjusted p-value
#res <- res[order(res$padj), ]
## Unir con los datos de counts
resdata <- (as.data.frame(res))#, as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
#names(resdata)[1] <- "Gene"
head(resdata)
## Escribir resultados
write.csv(resdata, file="normales_vs_stage1.csv")





