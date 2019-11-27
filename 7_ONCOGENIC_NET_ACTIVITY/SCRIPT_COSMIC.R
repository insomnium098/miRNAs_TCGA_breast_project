
library(dplyr)
library(matrixStats)
library(plotly)
library(gplots)

datos_census <- read.csv("Census_allWed Nov  6 20_29_51 2019.csv",
                         header = TRUE)
datos_network <- read.csv("DATASET_MIRNAS_SIGNIF_COOL.csv",
                          row.names = 1)


filtro <- filter(datos_network, Target %in% datos_census$Gene.Symbol)


###merge

union <- merge(filtro, datos_census, by.x = "Target",
               by.y = "Gene.Symbol")


######Hacer melt para hacer un heatmap

mirna <- "hsa-miR-106b-5p"

for (mirna in unique(as.character(union$MIRNA))){
  filt1 <- filter(union, MIRNA == mirna)
  #total_genes <- length(rownames(filt1))
  fin <- as.data.frame(mirna)
  ##fusion
  #filt_fusion <- filter(filt1, Role.in.Cancer == "fusion")
  #fusion <- length(rownames(filt_fusion))
  
  ######GENES QUE AMBOS SON ONCO Y TSG
  
  filt_ambos <- filter(filt1, Role.in.Cancer == "oncogene, TSG")
  ambos <- length(rownames(filt_ambos))
  fin$Oncogene_TSG <- ambos
  
  #fin$fusion <- fusion
  #####"oncogene"##### se toman en cuenta los fusion
  filt_oncogene_1 <- filter(filt1, Role.in.Cancer == "oncogene")
  
  filt_oncogene_2 <- filter(filt1, Role.in.Cancer == "oncogene, fusion")
  
  filt_oncogene <- rbind(filt_oncogene_1, filt_oncogene_2)
  
  ####INICIA BIEN V2
  ###Todos los oncogenes
  oncogene_all <- length(rownames(filt_oncogene))
  fin$Oncogene_all <-oncogene_all
  
  ###CASO 1: ONCOGENES RELACIONADOS POSITIVAMENTE
  ##TEST: FILTRAR SOLO AQUELLOS QUE TENGAN FC POSITIVO
  filt_oncogene_pos <- filt_oncogene[filt_oncogene$log2FoldChange >0,]
  oncogene_pos <- length(rownames(filt_oncogene_pos))
  
  fin$Oncogene_Up <- oncogene_pos
  ####CASO 2: ONCOGENES RELACIONADOS POSITIVAMENTE
  filt_oncogene_neg <- filt_oncogene[filt_oncogene$log2FoldChange <0,]
  oncogene_neg <- length(rownames(filt_oncogene_neg))
  fin$Oncogene_Down <- oncogene_neg
  
  
  ############TSG
  
  #"TSG"
  ##TEST: FILTRAR SOLO AQUELLOS QUE TENGAN FC negativo
  filt_tsg_1 <- filter(filt1, Role.in.Cancer == "TSG")
  filt_tsg_2 <- filter(filt1, Role.in.Cancer == "TSG, fusion")
  
  filt_tsg <- rbind(filt_tsg_1, filt_tsg_2)
  
  
  tsg_all <- length(rownames(filt_tsg))
  fin$TSG_all <- tsg_all
  
  ###CASO 1: TSG relacionado positivamente
  filt_tsg_pos<- filt_tsg[filt_tsg$log2FoldChange > 0,]
  tsg_pos <- length(rownames(filt_tsg_pos))
  fin$TSG_Up <- tsg_pos
  
  ###CASO 2: TSG RELACIONADO NEGATIVAMENTE
  
  filt_tsg_neg<- filt_tsg[filt_tsg$log2FoldChange < 0,]
  tsg_neg <- length(rownames(filt_tsg_neg))
  fin$TSG_Down <- tsg_neg
  ####TERMINA BIEN V2
  
  
  

  
  
  #########OLD DEPRECATED
  #"TSG"
  ##TEST: FILTRAR SOLO AQUELLOS QUE TENGAN FC negativo
  #filt_tsg <- filter(filt1, Role.in.Cancer == "TSG")
  #filt_tsg<- filt_tsg[filt_tsg$log2FoldChange <0,]
  
  #tsg <- length(rownames(filt_tsg))
  
  #fin$TSG <- tsg
  #"oncogene, fusion" 
  #filt_onco_fus <- filter(filt1, Role.in.Cancer == "oncogene, fusion")
  #onco_fus <- length(rownames(filt_onco_fus))
  
  #fin$Onco_Fus <- onco_fus
  #"oncogene, TSG"
  #filt_onco_tsg <- filter(filt1, Role.in.Cancer == "oncogene, TSG")
  #onco_tsg <- length(rownames(filt_onco_tsg))
  
  #fin$onco_tsg <- onco_tsg
  #"oncogene, TSG, fusion"
  #filt_onco_tsg_fus <- filter(filt1, Role.in.Cancer == "oncogene, TSG, fusion")
  #onco_tsg_fus <- length(rownames(filt_onco_tsg_fus))
  
  #fin$onco_tsg_fus <- onco_tsg_fus
  #"TSG, fusion" 
  #filt_tsg_fus <- filter(filt1, Role.in.Cancer == "TSG, fusion")
  #tsg_fus <- length(rownames(filt_tsg_fus))
  
  #fin$TSG_Fus <- tsg_fus
  
  
  total_genes <- oncogene_all + tsg_all + ambos#length(rownames(filt1))
  
  
  #Hallmark
  
  filt_hallmark <- filter(filt1, Hallmark == "Yes")
  hallmark <- length(rownames(filt_hallmark))
  
  fin$hallmark <- hallmark
  
  #MIRNA ONCOGENIC INDEX
  #V2: ONCO : onco_pos + tsg_neg
  #    TSG : onco_neg + tsg_pos
  ##es la resta de los onco menos tsg entre el total
  #
  #Valores cercanos a 1 significan oncogenicos y menores a 0 TSG
  ##A AMBOS INDEX SE LES RESTA LOS AMBOS
  
  index_onco <- (oncogene_pos + tsg_neg) - ambos
  index_tsg <- (oncogene_neg + tsg_pos) - ambos
  
  index <- (index_onco - index_tsg) / (total_genes)
  fin$index <- index
  
  
  ###REORDENAR DF PARA GRAFICAR MAS CHIDO
  
  #fin <- fin[,c(1,2,5,3,7,4,6,8,9)]
  
  fin <- fin[,c(1,2,3,6,4,8,5,7,9,10)]
  
  ###OLD
  ##es la resta de los onco menos tsg entre el total
  #
  #Valores cercanos a 1 significan oncogenicos y menores a 0 TSG
  ##No se tomaran en cuenta onco_tsg y onco_tsg_fus
  #index_onco <- oncogene + onco_fus
  #index_tsg <- tsg + tsg_fus
  
  #index <- (index_onco - index_tsg) / (total_genes)
  
  #fin$index <- index
  
  ###END OLD
  
  
  #Not annotated
  #filt_na
  #filt_na <- filter(filt1, Role.in.Cancer == "")
  #not_ann <- length(rownames(filt_na))
  
  #fin$not_ann <- not_ann
  #Si la nueva tabla no existe, se crea
  if (!exists("df_final")){
    df_final <- fin
  } else{
    #Si la nueva tabla existe, se une 
    temp_dataset <- fin
    df_final <- rbind(df_final,temp_dataset)
    rm(temp_dataset)
    
  }
  
  
  

  
}


#####REEMPLAZAR NAN POR 0

df_final[df_final$index == "NaN",9] <- 0

write.csv(df_final,"DF_FINAL_COSMIC.csv")



