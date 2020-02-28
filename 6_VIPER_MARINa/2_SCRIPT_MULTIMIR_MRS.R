
library(multiMiR)
library(dplyr)
library(plyr)


datos_network <- read.delim("network_SIN_DPI.txt",header = TRUE)

datos_mrs <- read.delim("mirnas.interesantes.txt", header = FALSE)

mirnas_todos <- as.character(unique(datos_mrs$V1)) 


if(exists("dataset_mirnas")){
  mirnas_todos <- mirnas_todos[counter:length(mirnas_todos)]
  
} else{
  counter <- 0
}

#mirna_nombre <- "hsa-miR-18a-5p"


for (mirna_nombre in mirnas_todos){
  print(paste0("Verificando: ", mirna_nombre))
  
  counter <- counter + 1
  cat("\r", counter, "de", length(mirnas_todos),"miRNAs", "\r")
  
  filtro_mirna <- filter(datos_network, Regulator == mirna_nombre)
  genes_mirna_filtro <- as.character(filtro_mirna$Target)
  
  data_multimir <- get_multimir(org     = "hsa",
                                mirna  = mirna_nombre,
                                target = genes_mirna_filtro,
                                table   = c("all"),
                                summary = FALSE,
                                predicted.cutoff      = 300000,
                                predicted.cutoff.type = "n",
                                predicted.site        = "all")
  
  resultados <- (data_multimir@data)
  ###AQUI ERROR
  resultados <- filter(resultados, type == "predicted" |
                         type == "validated")
  

  
  

  
  
  
  if (length(rownames(resultados))== 0){
    if (!exists("dataset_mirnas_sin_blancos")){
      dataset_mirnas_sin_blancos <- mirna_nombre
    } else {
      
      temp_dataset <- mirna_nombre
      dataset_mirnas_sin_blancos <- c(dataset_mirnas_sin_blancos,
                                      temp_dataset)
      rm(temp_dataset)
      paste0(mirna_nombre," no tiene blancos")
      
    }
    
    next
  } else{
    
    
    
    #########
    
    ###Filtrar los predichos de netwrok para
    ### extraer mi
    
    
    
    filtro_mi <- filter(filtro_mirna, Target %in% resultados$target_symbol)
    
    ###Extraer unicos de resultados
    
    
    resultados_unique <- resultados[!duplicated(
      resultados$target_symbol
    ),]
    
    
    ###Remover genes de resultados_unique que no esten en 
    ### filtro_mi
    
    resultados_unique <- filter(resultados_unique, target_symbol %in%
                                  filtro_mi$Target)
    
    ####Unir los dos ultimos df
    ####AQUI ERROR, HACER MERGE
    
    #resultados <- cbind(filtro_mi, resultados_unique)
    resultados <- merge(filtro_mi, resultados_unique,
                  by.x = "Target", by.y = "target_symbol")
    #######
    #Si el datasset no existe, se crea
    if (!exists("dataset_mirnas")){
      dataset_mirnas <- resultados
    } else {
      
      temp_dataset <- resultados
      
      ###
      dataset_mirnas <- rbind.fill(dataset_mirnas,temp_dataset)
      rm(temp_dataset)
      
    }
  }
  
}


write.csv(dataset_mirnas,"PREDICCION_y_VALIDACION_MRS_SIN_SPI_CORREGIDO.csv")

#write.table(dataset_mirnas[,1:3],"network_predicted_SIN_DPI.txt")





