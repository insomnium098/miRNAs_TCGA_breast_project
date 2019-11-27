

library(dplyr)
######SOLO SE HARA LA CORRELACION EN TUMORES


####
tp_genes <- read.csv("diffexpr-results_DESEQ2.csv",header = TRUE, row.names = 1)
  
####EXTRAER SOLO SIGNIFICATIVOS

tp_genes_filter <- filter(tp_genes, padj < 0.01)

tp_genes <- tp_genes_filter
nombres_nt <- read.csv("matriz_NT_GENES.csv",header = TRUE)


extraccion_tumores <- tp_genes[,!colnames(tp_genes) %in% colnames(nombres_nt)]
rownames(extraccion_tumores) <- make.unique(as.character(extraccion_tumores$Gene))
extraccion_tumores <- extraccion_tumores[,8:length(colnames(extraccion_tumores))]

tp_genes <- extraccion_tumores
#tp_genes <- tp_genes[,100:length(colnames(tp_genes))]
#nombres_tp <- as.character(colnames(tp_genes))
#nombres_tp <- gsub(".1","",nombres_tp, fixed = TRUE)
####
##TOMAR EN CUENTA SOLO LOS PRIMEROS 12 CARACTERES DE LOS COLNAMES
#nombres_tp <- substr(nombres_tp, 0, 12)
#colnames(tp_genes) <- nombres_tp

############EXTRAER SOLO LOS SIGNIFICATIVOS

#diff_exp_archivo <- read.csv("RESULTADOS_SUBEXP_STAGE1.csv",header = TRUE,
#                             row.names = 1)


#nombres_extraccion_diff_exp <- as.character(diff_exp_archivo$X)

#tp_genes <- rownames(tp_genes)

####################
tp_mirnas <- read.csv("vst_transformation_mirnas.csv", header = TRUE, row.names = 1)
tp_mirnas <- tp_mirnas[,90:length(colnames(tp_mirnas))]

nombres_tp_mirnas <- as.character(colnames(tp_mirnas))
nombres_tp_mirnas <- gsub("stage1.","",nombres_tp_mirnas, fixed = TRUE)
nombres_tp_mirnas <- gsub("stage2.","",nombres_tp_mirnas, fixed = TRUE)
nombres_tp_mirnas <- gsub("stage3.","",nombres_tp_mirnas, fixed = TRUE)
nombres_tp_mirnas <- gsub("stage4.","",nombres_tp_mirnas, fixed = TRUE)



colnames(tp_mirnas) <- nombres_tp_mirnas

mirnas_tp <- as.character(rownames(tp_mirnas))
mirnas_tp <- gsub(".","-",mirnas_tp,fixed = TRUE)
rownames(tp_mirnas) <- mirnas_tp

 #####EXTRAER SOLO MIRNAS de los resultados

resultados_mirnas <- read.csv("RESULTADOS_MIRNAS_SOBRE_FC_0_5.csv", header = TRUE,
                              row.names = 1)

mirnas_interesantes <- as.character(resultados_mirnas$X)

tp_mirnas <- tp_mirnas[mirnas_interesantes,]

##############


nombres_tp_genes <- colnames(tp_genes)

nombres_tp_mirnas <- colnames(tp_mirnas)


########INTERSECCION TUMORES

nombres_tumores_intersect <- intersect(nombres_tp_genes,nombres_tp_mirnas)

############EXTRACCION



##TUMORES
matriz_tumores_genes_intersect <- tp_genes[,nombres_tumores_intersect]
matriz_tumores_mirnas_intersect <- tp_mirnas[,nombres_tumores_intersect]
matriz_tumores_intersect <- rbind(matriz_tumores_genes_intersect,matriz_tumores_mirnas_intersect)


genes_tumores_intersect <- rownames(matriz_tumores_intersect)

library(dplyr)

matriz_tumores_intersect <- mutate_all(matriz_tumores_intersect, function(x) as.numeric(as.character(x)))

rownames(matriz_tumores_intersect) <- genes_tumores_intersect

matriz_tumores_intersect_2 <- cbind(genes_tumores_intersect,matriz_tumores_intersect)
colnames(matriz_tumores_intersect_2)[1] <- "gene"
#####UNION FINAL

#rm(list=ls(all=TRUE))
write.table(matriz_tumores_intersect_2,"MATRIZ_UNION_FINAL_MIT_ARACNE.txt"
            ,sep = "\t",row.names = FALSE,quote = FALSE)

write.table(mirnas_interesantes,"mirnas.interesantes.txt",
            row.names = FALSE,col.names = FALSE,
            quote = FALSE)



