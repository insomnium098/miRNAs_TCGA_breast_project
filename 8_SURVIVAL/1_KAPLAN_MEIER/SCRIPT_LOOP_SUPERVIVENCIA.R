rm(list = ls(all.names = TRUE)) 
library(survival)
library(survminer)
library(dplyr)
####NOTA#####
######EL SCRIPT SOLO FUNCIONA CON UN ESTADIO CLINICO A LA VEZ
#######
#######VERSION 2.0, SE AGREGO EL FILTRO DE DATOS CLINICOS
########





datos_clinicos <- read.csv("COMPLETOS_DATOS_CLINICOS.csv", header = TRUE, row.names = 2)
#############FILTRO DE DATOS CLINICOS

###Censura a 5 años
#5 AÑOS

##ALL YEARS 9125
dias_filtro <- c(9125)

if(exists("dias_filtro")==TRUE){
  print("EXISTE CENSURA DE DIAS")
} else{
  print("NO EXISTE CENSURA DE DIAS")
  #dias_filtro <- max(days_followup)
}

days_followup <- datos_clinicos$days_to_last_follow_up

#days_followup[days_followup > 1825] <- 1825
days_followup <- days_followup / 365

days_to_death <- datos_clinicos$days_to_death
#days_to_death[days_to_death > 1825] <- 1825
days_to_death<- days_to_death / 365

datos_clinicos$days_to_last_follow_up <- days_followup
datos_clinicos$days_to_death <- days_to_death

#datos_clinicos_t <- filter(datos_clinicos, days_to_last_follow_up <= dias_filtro)
#rownames(datos_clinicos_t) <- as.character(datos_clinicos_t$submitter_id)

#datos_clinicos <- datos_clinicos_t

####NOTA#####
######EL SCRIPT SOLO FUNCIONA CON UN ESTADIO CLINICO A LA VEZ
#######


#######################

###CONFIGURACION_LOOP
#estadios <- c(4)

datos_mirnas <- read.csv("MRS_SIGNIF_05.csv",
                         header = TRUE)

mirnas_loop <- as.character(datos_mirnas$Regulon)

#mirnas_loop <- as.character(mirnas_cool$X)
#mirnas_loop <- gsub(".","-", mirnas_loop, fixed = TRUE)

##CUANTILES EXPRESION
threshold_up <-0.50

##SI SE USA EL THRESHHOLD DOWN MODIFICAR LA COMPARACION DE CUANTILES EN EL LOOP
#threshold_down <- 0.33
#####################

nombre_archivo_estadio <- c("vst_transformation.csv")#read.csv("mariz_tp_brca_readcount.csv", header = TRUE, row.names = 1)
#######Remover a los pacientes sanos


###en esta ,atriz estan todos con rld

rld <- read.csv("vst_transformation.csv", header = TRUE, row.names = 1,
                stringsAsFactors = FALSE)

#######Remover a los pacientes sanos
rld <- rld[,-c(1:90)]

nombres_mal <- colnames(rld)
nombres_bien <- gsub(".","-", nombres_mal, fixed = TRUE)
nombres_bien <- gsub("stage1-","", nombres_bien, fixed = TRUE)
nombres_bien <- gsub("stage2-","", nombres_bien, fixed = TRUE)
nombres_bien <- gsub("stage3-","", nombres_bien, fixed = TRUE)
nombres_bien <- gsub("stage4-","", nombres_bien, fixed = TRUE)

colnames(rld) <- nombres_bien

rld_t <- rld#as.data.frame(t(rld))
#################### INICIO CONFIGURACION
#for (estadio_clinico in estadios){
  #print (estadio_clinico)

counter <- 0
    for (microrna in mirnas_loop) {
    mirna <- microrna#c("hsa-mir-592")
    #threshold_up <-0.50
    #threshold_down <- 0.33
    counter <- counter + 1
    
    longitud_mirnas <- length(mirnas_loop)
    print(paste0(c(counter," de ",longitud_mirnas)))
    ##########FIN DE CONFIGURACION
    
    #######Extraer mirna de RLD
    #BIEN
    rld_extraido <- rld_t[mirna,]
    rld_extraido <- as.data.frame(sapply(rld_extraido , as.numeric))
    colnames(rld_extraido) <- mirna
    rld_extraido <- as.data.frame(t(rld_extraido))
    colnames(rld_extraido) <- colnames(rld)
    
    ######extraer stage de rld con mirna
    
    #tcga_stage_nombres <- paste0("nombres_stage",estadio_clinico)
    #print(tcga_stage_nombres)
    
    rld_mirna_stage <- rld_extraido#[, get(tcga_stage_nombres)]
    
    ##########voltear extraccion
    #BIEN
    rld_mirna_stage_t <- as.data.frame(t(rld_mirna_stage))
    rownames(rld_mirna_stage_t) <- gsub(".","-",rownames(rld_mirna_stage_t),
                                        fixed = TRUE)
    
    ##########
    #########EXTRAER SOLO AQUELLOS PACIENTES CON DATOS CLINICOS
    pacientes_con_clinicos <- as.character(rownames(datos_clinicos))
    pacientes_mirnaseq_mal <- colnames(rld)
    ###remover duplicados con string -1
    #pacientes_duplicados <- grep("-1",pacientes_mirnaseq_mal)
    #pacientes_duplicados_2 <- grep("-2",pacientes_mirnaseq_mal)
    #eliminacion_duplicados <- gsub("-1","",pacientes_mirnaseq_mal,fixed = TRUE)
    #nombres_mirnaseq_bien <- pacientes_mirnaseq_mal[-pacientes_duplicados]
    nombres_mirnaseq_bien_index <- pacientes_mirnaseq_mal %in%pacientes_con_clinicos#unique(pacientes_con_clinicos)#nombres_mirnaseq_bien[-pacientes_duplicados_2]
    nombres_mirnaseq_bien <- pacientes_mirnaseq_mal[nombres_mirnaseq_bien_index]
  
    
    
    
    
    ################FUNCION CALCULO DE EXPRESION EN CUARTILES
    if (exists("threshold_down")){
    cuantiles <- quantile(rld_mirna_stage_t[,mirna], c(threshold_up, threshold_down))
    valor_cuantil_up <- as.numeric(cuantiles[1])
    valor_cuantil_down <- as.numeric(cuantiles[2])
    }
    else{
    cuantiles <- quantile(rld_mirna_stage_t[,mirna], c(threshold_up))
    valor_cuantil_up <- as.numeric(cuantiles[1])
    }
    
    cuantil_normal <- quantile(rld_mirna_stage_t[,mirna], c(0.60))
    
    
    ####En algunos mirnas, como mir-2114-3p la mediana es el 
    #mismo valor que el mas bajo, en esos casos se usara el promedio
    
    valor_menor <- min(rld_mirna_stage_t[,mirna])
    
    valor_cuantil_up <- as.numeric(cuantiles[1])
    
    if(valor_menor == valor_cuantil_up){
      valor_cuantil_up <- mean(rld_mirna_stage_t[,mirna])
      if (!exists("dataset_mirnas_con_mean")){
        dataset_mirnas_con_mean <- as.data.frame(mirna)
      } else{
        #Si el dataset existe, se une 
        temp_dataset <- as.data.frame(mirna)
        dataset_mirnas_con_mean <- rbind(dataset_mirnas_con_mean,temp_dataset)
        rm(temp_dataset)
        
      }
      
    }
    
    valor_cuantil_normal <- as.numeric(cuantil_normal[1])
    
    ############FUNCION EXTRACCION EN BASE A CUANTILES
    rld_mirna_stage_t$grupo_expresion <- c()

    #rld_mirna_stage_t[rld_mirna_stage_t$grupo_expresion > valor_cuantil_up] <- 5
    for (paciente in (nombres_mirnaseq_bien)){
      #for (valor in (rld_mirna_stage_t[,mirna])) {
      valor <- rld_mirna_stage_t[paciente,]
      if (valor >= valor_cuantil_up){
        expresion <- c("SOBREXPRESADO")
        #} else if (valor == valor_cuantil_normal){
        #  expresion <- c("EXPRESION_NORMAL")
      } else if (valor < valor_cuantil_up){
        expresion <- c("SUBEXPRESADO")
      } 
      
      #####añadir las diferentes columnas
      expresion <- as.data.frame(expresion)
      rownames(expresion) <- paciente
      expresion$valor_expresion <- valor
      ###añadir mirna
      expresion$microrna <- mirna
      #expresion <- c(valor,expresion)
      #Si el datasset no existe, se crea
      if (!exists("dataset_expresion")){
        dataset_expresion <- expresion 
      }
      
      
      #Si el dataset existe, se une 
      if (exists("dataset_expresion")){
        temp_dataset <-expresion
        #dataset<-rbind(dataset, temp_dataset)
        dataset_expresion <- rbind(dataset_expresion,temp_dataset)
        rm(temp_dataset)
      }
      
      
      ####extraer solo los pacientes originales
      nombres_originales <- nombres_mirnaseq_bien#colnames(rld_mirna_stage)
      ###remover duplicados con string -1
      #pacientes_duplicados <- grep("-1",nombres_originales)
      #nombres_originales_bien <- nombres_originales[-pacientes_duplicados]
      dataset_final <- dataset_expresion[nombres_originales,]
    }
    
    #con este comando se muestra cual es el paciente que se repite por el loop
    #dataset_mal <- dataset_expresion[!(rownames(dataset_expresion) %in% nombres_originales), ]
    
    #######AÑADIR DATOS CLINICOS
    
    #library(TCGAbiolinks)
    # Survival Analysis SA
    #clinical_patient_Cancer <- GDCquery_clinic("TCGA-BRCA","clinical")
    
    ###LOOP PARA EXTRAER EL VITAL STATUS Y REEMPLAZAR ALIVE POR 0 Y DEAD POR 1
    
    for (paciente in (nombres_mirnaseq_bien)){
      #for (valor in (rld_mirna_stage_t[,mirna])) {
      status_vital <- datos_clinicos[paciente,"vital_status"]
      if (is.na(status_vital)){
        break
      } else if (status_vital == c("dead")){
        status_vital_numero <- c(1)
      } else if (status_vital == c("alive")){
        status_vital_numero <- c(0)
      }
      
      #####añadir las diferentes columnas
      status_vital_numero  <- as.data.frame(status_vital_numero)
      rownames(status_vital_numero) <- paciente
      
      #Si el datasset no existe, se crea
      if (!exists("dataset_status_vital")){
        dataset_status_vital <- status_vital_numero 
      }
      
      
      #Si el dataset existe, se une 
      if (exists("dataset_status_vital")){
        temp_dataset <-status_vital_numero
        #dataset<-rbind(dataset, temp_dataset)
        dataset_status_vital <- rbind(dataset_status_vital,temp_dataset)
        rm(temp_dataset)
      }
      
      
      
      ####extraer solo los pacientes originales
      nombres_originales <- nombres_mirnaseq_bien#colnames(rld_mirna_stage)
      dataset_status_vital_final <- as.data.frame(dataset_status_vital[nombres_originales,])
      rownames(dataset_status_vital_final) <- nombres_originales
      colnames(dataset_status_vital_final) <- c("status_vital")
    }

  ######EXTRAER MAS DATOS CLINICOS
  #days to death
  dias_a_muerte <- as.data.frame(datos_clinicos$days_to_death)
  rownames(dias_a_muerte) <- rownames(datos_clinicos)
  dias_a_muerte <- as.data.frame(dias_a_muerte[nombres_mirnaseq_bien,])
  rownames(dias_a_muerte) <- nombres_mirnaseq_bien
  colnames(dias_a_muerte) <- c("days_to_death")
  dias_a_muerte[is.na(dias_a_muerte)] <- 0
  ##
  
  #days to last followup
  dias_followup <- as.data.frame(datos_clinicos$days_to_last_follow_up)
  rownames(dias_followup) <- rownames(datos_clinicos)
  dias_followup <- as.data.frame(dias_followup[nombres_mirnaseq_bien,])
  rownames(dias_followup) <- nombres_mirnaseq_bien
  colnames(dias_followup) <- c("days_to_last_followup")
  dias_followup[is.na(dias_followup)] <- 0
  ##
  
  
  
  ######UNIR TODAS LAS MATRICES
  
  dataset_chingon <- cbind(dataset_final, dataset_status_vital_final, dias_a_muerte, dias_followup)

########CREAR COLUMNA DIAS, SI EL PACIENTE MURIO SE USA EL DAYS TO DEATH, EN CASO CONTRARIO SE USA EL
#######LAST FOLLOWUP

for (paciente in (nombres_mirnaseq_bien)){
  #for (valor in (rld_mirna_stage_t[,mirna])) {
  status_vital <- dataset_chingon[paciente,"status_vital"]
  if (status_vital == c("1")){
    dias_bien <- dataset_chingon[paciente,"days_to_death"]
  } else if (status_vital == c("0")){
    dias_bien <- dataset_chingon[paciente,"days_to_last_followup"]
  } 
  
  #####añadir las diferentes columnas
  dias_bien_df  <- as.data.frame(dias_bien)
  rownames(dias_bien_df) <- paciente
  
  #Si el datasset no existe, se crea
  if (!exists("dataset_dias_bien")){
    dataset_dias_bien <- dias_bien_df 
  }
  
  
  #Si el dataset existe, se une 
  if (exists("dataset_dias_bien")){
    temp_dataset <-dias_bien_df
    #dataset<-rbind(dataset, temp_dataset)
    dataset_dias_bien <- rbind(dataset_dias_bien,temp_dataset)
    rm(temp_dataset)
  }
  
  ####extraer solo los pacientes originales
  nombres_originales <- nombres_mirnaseq_bien#as.character(rownames(datos_clinicos))
  dataset_dias_bien_final <- as.data.frame(dataset_dias_bien[nombres_originales,])
  rownames(dataset_dias_bien_final) <- nombres_originales
  colnames(dataset_dias_bien_final) <- c("dias_bien")
}

dataset_Super_chingon <- cbind(dataset_chingon, dataset_dias_bien_final)

#correcion_estadio_clinico <- gsub("nombres_","",tcga_stage_nombres)
nombre_archivo_final <- paste0("datos_clinicos","_","_",microrna,".csv")

#titulo_archivo_chingon <- paste0(mirna,"-","stage-",estadio_clinico,".csv")

####ESCRIBIR DATOS CLINICOS CON EXPRESION
#write.csv(dataset_Super_chingon, nombre_archivo_final)



#####OBTENER PVALUES DE KAPLAN MEIER

log_rank <- survdiff(Surv(dias_bien, status_vital) ~ expresion, data = dataset_Super_chingon)
p.val <- 1 - pchisq(log_rank$chisq, length(log_rank$n) - 1)

#Si el datasset no existe, se crea
if (!exists("dataset_p_values")){
  dataset_p_values <- as.data.frame(p.val)
  dataset_p_values$microrna <- mirna
} else{
  #Si el dataset existe, se une 
  temp_dataset <- as.data.frame(p.val)
  temp_dataset$microrna <- mirna
  dataset_p_values <- rbind(dataset_p_values,temp_dataset)
  rm(temp_dataset)
  
}

#####Si es significativo entonces se grafica

if (p.val < 0.05){

#AQUI PONER CODIGO DE KAPLAN MEIER
###TITULO
titulo_grafico <- paste0(mirna)
titulo_archivo <- paste0(mirna,".pdf")

#######RENIVELAR LOS FACTORES

dataset_Super_chingon$expresion <- relevel(dataset_Super_chingon$expresion, ref = "SUBEXPRESADO")
###########

fit_bien <- survfit(Surv(dias_bien, status_vital) ~ expresion, data = dataset_Super_chingon)
#plota <- ggsurvplot(fit, palette = "grey",legend.labs = levels(dataset_Super_chingon$expresion),
#                    pval = TRUE, data = dataset_Super_chingon, censor.shape="|", censor.size = 20, title=titulo_grafico,
#                    pval.method = TRUE, xlim = c(0,dias_filtro), ggtheme = theme_minimal()
#                    )

plota <- ggsurvplot(fit = fit_bien, palette = "grey",#c("black", "black"),
                    pval = TRUE, censor.shape="|", censor.size = 3, #title=titulo_grafico,
                    pval.method = FALSE, data = dataset_Super_chingon,
                    linetype = "strata",
                    #legend = c(0.8,0.3), 
                    legend.labs = c("Low", "High"),
                    legend.title = mirna,
                    #xlab = "Time (months)",
                    font.x = "bold",
                    title= "TCGA all stages",
                    font.main = "bold",
                    font.submain = "bold",
                    font.caption = "bold",
                    font.y = "bold",
                    font.tickslab = "bold",
                    font.legend = "bold",
                    surv.scale = c("default"),
                    ylab="Overall survival (%)",
                    xlab="Time (Years)",
                    fun = function(y) y*100,
                    pval.coord = c(3, 21),
                    pval.method.coord = c(3, 28),
                    xlim = c(1,(dias_filtro/365)))

###Centrar titulo
plota$plot <- plota$plot +  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

######GENERAR imagen
pdf(titulo_archivo, onefile = FALSE)
print(plota)
dev.off()

}

##############

########ELIMINAR DATASET PARA CONTINUAR CON EL LOOP
rm(dataset_chingon,dataset_dias_bien,dataset_dias_bien_final,dataset_expresion,
   dataset_final,dataset_status_vital,dataset_status_vital_final,
   dataset_Super_chingon)
################
}
  #print(estadios)###MAL
  #print(nombre_archivo_final)
  #print(nrow(dataset_Super_chingon))
#}

write.csv(dataset_p_values,"DATOS_P_VALUES_KAPLAN.csv")
write.csv(dataset_mirnas_con_mean,"mirnas_con_mean.csv")


