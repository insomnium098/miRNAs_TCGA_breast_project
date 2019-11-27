
library(dplyr)
library(org.Hs.eg.db)
library(annotate)

###DATOS DE MRNA EN DESEQ2 RLOG
#mrna_deseq2 <-  "/Users/daniel/OneDrive\ -\ UNIVERSIDAD\ NACIONAL\ AUTÓNOMA\ DE\ MÉXICO/owncloud/X99/TCGA-BRCA/DIFF\ EXP\ MRNA\ BRCA\ NORMALES\ VS\ TUMORALES/diffexpr-results_DESEQ2.csv" 
#datos_mrna_deseq2 <- read.csv(mrna_deseq2, header = TRUE)

#PACIENTES_NT MRNA PARA NAMES
datos_mrna_nt <- read.csv("../../DATOS_PROTEIN_CODING/MATRIZ_VST_PROTEIN_CODING_NT.csv", header = TRUE,
                          row.names = 1)
datos_mrna_nt_exp <- datos_mrna_nt
datos_mrna_nt <- colnames(datos_mrna_nt)#[2:length(colnames(datos_mrna_nt))]

#PACIENTES_TP MRNA PARA NAMES
datos_mrna_tp <- read.csv("../../DATOS_PROTEIN_CODING/MATRIZ_VST_PROTEIN_CODING_TP.csv", header = TRUE,
                          row.names = 1)
datos_mrna_tp_exp <- datos_mrna_tp
datos_mrna_tp <- colnames(datos_mrna_tp)#[2:length(colnames(datos_mrna_tp))]

####Filtro_signif

#filtro_signif <-datos_mrna_deseq2 #filter(datos_mrna_deseq2, padj < 0.05)


#######EXTRAER SOLO AQUELLOS GENES PRESENTES Y SIGNIFICATIVOS EN LA RED

#aracne_network <- read.delim("network.txt", header = TRUE)
#genes_aracne <- as.character(unique(aracne_network$Target))

#filtro_aracne <- filter(filtro_signif, Gene %in% genes_aracne)

#filtro_signif <- filtro_aracne

#############DATOS DE MIRNAS
mirnas_chidos <- read.csv("../../DATOS_ANALISIS/RESULTADOS_MIRNAS_SOBRE_FC_0_5.csv",
                          header = TRUE, row.names = 1)
mirnas_chidos <- as.character(mirnas_chidos$X)
mirnas_vst <- read.csv("../../DATOS_ANALISIS/vst_transformation_mirnas.csv",header = TRUE,
                       row.names = 1)

nombres_mirnas_vst_mal <- colnames(mirnas_vst)
nombres_mirnas_vst_mal <- gsub("stage1.","",nombres_mirnas_vst_mal)
nombres_mirnas_vst_mal <- gsub("stage2.","",nombres_mirnas_vst_mal)
nombres_mirnas_vst_mal <- gsub("stage3.","",nombres_mirnas_vst_mal)
nombres_mirnas_vst_mal <- gsub("stage4.","",nombres_mirnas_vst_mal)
colnames(mirnas_vst) <- nombres_mirnas_vst_mal




mirnas_nt <- mirnas_vst[,1:90]
mirnas_tp <- mirnas_vst[,91:length(colnames(mirnas_vst))]

mirnas_vst_extraccion <- mirnas_vst[rownames(mirnas_vst) 
                                             %in% mirnas_chidos,]

#mirnas_vst_normales <- mirnas_vst_extraccion[,1:90]

#mirnas_vst_tp <- mirnas_vst_extraccion[,90:length(colnames(mirnas_vst_extraccion))]

######INTERSECT DE PACIENTES
###NT

pacientes_nt_intersect <- intersect(datos_mrna_nt,
                                    colnames(mirnas_nt))

pacientes_tp_intersect <- intersect(datos_mrna_tp,
                                    colnames(mirnas_tp))

#####EXTRACCION DE PACIENTES DE MRNA

#rownames(filtro_signif) <- as.character(filtro_signif$Gene)

nt_mrna_intersect <- datos_mrna_nt_exp[,colnames(datos_mrna_nt_exp) %in%
                                     pacientes_nt_intersect]

tp_mrna_intersect <- datos_mrna_tp_exp[,colnames(datos_mrna_tp_exp) %in%
                                     pacientes_tp_intersect]

######EXTRACCION DE PACIENTES DE MIRNAS

nt_mirnas_intersect <- mirnas_vst_extraccion[,
                                             colnames(mirnas_vst_extraccion)[1:90] %in% pacientes_nt_intersect]

nt_mirnas_intersect <- nt_mirnas_intersect[,1:90]



vst_tp <- mirnas_vst_extraccion[,91:length(colnames(mirnas_vst_extraccion))]
tp_mirnas_intersect <- vst_tp[,colnames(vst_tp)
                              %in% pacientes_tp_intersect]

######Convertir ENSEMBL A GENESYMBOL 
conversor_id_ensemble <- function(matriz)
{
  datos_conversor <- matriz #c("BC069070")#read.delim("GENESIS_TABLA_SIGNIFICATIVOS.txt",header= TRUE, row.names = 1)
  resultados_symbol = mapIds(org.Hs.eg.db,
                             keys=rownames(datos_conversor), 
                             column="SYMBOL",#En esta opcion se selecciona que valor se quiere obtener de la conversion
                             keytype="ENSEMBL",#En esta opcion se selecciona en que formato estan los datos de origen
                             multiVals="first")
  #Reportar resultados
  converios_data_frame <- as.data.frame(resultados_symbol)
  conversion_bien_bind <- cbind(converios_data_frame, matriz)
  Resultado_final <- conversion_bien_bind#[complete.cases(conversion_bien_bind),]
  #row.names(matriz) <- resultados_symbol
  return(Resultado_final)
}

#################END CONVERSION IDS


conversion_tp <- conversor_id_ensemble(tp_mrna_intersect)
conversion_nt <- conversor_id_ensemble(nt_mrna_intersect)

conversion_tp <- conversion_tp[!is.na(conversion_tp$resultados_symbol),]
conversion_nt <- conversion_nt[!is.na(conversion_nt$resultados_symbol),]


rownames(conversion_tp) <- make.unique(as.character(conversion_tp$resultados_symbol))
rownames(conversion_nt) <- make.unique(as.character(conversion_nt$resultados_symbol))

conversion_tp <- conversion_tp[,2:length(colnames(conversion_tp))]
conversion_nt <- conversion_nt[,2:length(colnames(conversion_nt))]

tp_mrna_intersect <- conversion_tp
nt_mrna_intersect <- conversion_nt


##########UNION_NT_FINAL

matriz_union_nt_final <- rbind(nt_mrna_intersect,
                               nt_mirnas_intersect)


matriz_union_tp_final <- rbind(tp_mrna_intersect,
                               tp_mirnas_intersect)

matriz_union_nt_final$NAMES <- rownames(matriz_union_nt_final)
matriz_union_tp_final$NAMES <- rownames(matriz_union_tp_final)
#####Eliminar aquellos que no esten en genes_Aracne


nt_bien <- matriz_union_nt_final#filter(matriz_union_nt_final, 
          #        rownames(matriz_union_nt_final) %in% genes_aracne)
tp_bien <- matriz_union_tp_final#filter(matriz_union_tp_final, 
            #      rownames(matriz_union_tp_final) %in% genes_aracne)

rownames(nt_bien) <- nt_bien$NAMES
rownames(tp_bien) <- tp_bien$NAMES

nt_bien <- subset(nt_bien,select=-c(NAMES))
tp_bien <- subset(tp_bien,select=-c(NAMES))

write.csv(nt_bien,"MATRIZ_NT_VIPER_FINAL.csv")
write.csv(tp_bien,"MATRIZ_TP_VIPER_FINAL.csv")

