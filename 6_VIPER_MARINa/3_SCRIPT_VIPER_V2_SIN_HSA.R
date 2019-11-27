library(bcellViper)
library(mixtools)
library(viper)
library(SummarizedExperiment)
library(dplyr)


#data(bcellViper)



network <- read.delim("../network_predicted_SIN_DPI.txt",header = TRUE)
network <- network[,1:3]

###REMOVER HSA
network$Regulator <- gsub("hsa-","",
                          as.character(network$Regulator),
                          fixed = TRUE)

##FIX SEP AND AGO GENES

sep_genes <- network
sep_genes$Target <- gsub("sep-","SEPT", as.character(sep_genes$Target),
                         fixed = TRUE)
sep_genes$Target <- gsub("SEPT0","SEPT", as.character(sep_genes$Target),
                         fixed = TRUE)
sep_genes$Target <- gsub("ago-01","AGO1", as.character(sep_genes$Target),
                         fixed = TRUE)
network <- sep_genes

#####ELIMINAR HEADER DEL ARCHIVO
write.table(network,"network_viper.txt",
            sep = "\t",row.names = FALSE,quote = FALSE,
            col.names = FALSE)

adjfile <- file.path("network_viper.txt")
#####Datos de expresion con pacientes sanos y tumorales
##### normalizados con vst
###Archivo de diff_exp
matriz_nt <- read.csv("../MATRIZ_NT_VIPER_FINAL.csv",header = TRUE, row.names = 1)
matriz_tp <- read.csv("../MATRIZ_TP_VIPER_FINAL.csv",header = TRUE, row.names = 1)
### REMOVER HSA
row_nt <- rownames(matriz_nt)
row_nt <- gsub("hsa-","",row_nt,fixed = TRUE)
rownames(matriz_nt) <- row_nt

row_tp <- rownames(matriz_tp)
row_tp <- gsub("hsa-","",row_tp,fixed = TRUE)
rownames(matriz_tp) <- row_tp


union_matrices <- as.matrix(cbind(matriz_nt,matriz_tp))

#Datos_diff_exp <- file.choose()



regul <- viper::aracne2regulon(adjfile, union_matrices, format = "3col",
                               verbose = TRUE)




norm <- as.matrix(matriz_nt)
test <- as.matrix(matriz_tp)

signature <- rowTtest(norm, test)

signature_2 <- (qnorm(signature$p.value/2, lower.tail=F) * sign(signature$statistic))[, 1]

if (file.exists("NULL_MODEL.rds") == TRUE){
    nullmodel <- readRDS("NULL_MODEL.rds")
} else{
  nullmodel <- viper::ttestNull(norm, test, per=1000, repos=T, cores = 7,
                                verbose = TRUE)
  saveRDS(nullmodel,"NULL_MODEL.rds")
  
}

###Quitar HSA DE NULL MODEL

rows_null <- rownames(nullmodel)
rows_null <- gsub("hsa-","",rows_null,fixed = TRUE)
rownames(nullmodel) <- rows_null

######V1 CHIDA CON MARINA
#library(ssmarina)

#mrs <- marina(signature_2, regul, nullmodel)
#summary(mrs)

#ledger <- ledge(mrs)
#summary(ledger)


#mrshadow <- shadow(mrs, pval=25)
#combinatorial <- marinaCombinatorial(mrs, regulators=25)

######V2 CON VIPER
mrs_2 <- msviper(signature_2, regul, nullmodel, cores = 7,
                 verbose= TRUE)

test <- summary(mrs_2,217)
test <- test[with(test, order(FDR, -Size)), ]

plot(mrs_2, 50,cex = 1)

mrs_ledge_2 <- viper::ledge(mrs_2)
summary(mrs_ledge_2)

test <- summary(mrs_ledge_2,217)
##obtener ledges de cada miRNA
test2 <- mrs_ledge_2$ledge

test <- summary(mrs_ledge_2,217)
test <- test[with(test, order(FDR, -Size)), ]

####### Extraer los miRNAS con pval < 0.05
###### Para predecir sus blancos con mutimir

###

mirnas_signif <- as.character(mrs_signif$Regulon)


######
ledge_extraccion <- mrs_ledge_2$regulon

for (mirna_significativo in mirnas_signif){
  print(mirna_significativo)
  ledge_extraccion_1 <- ledge_extraccion[mirna_significativo]
  targets_ledge <- names(ledge_extraccion_1[[mirna_significativo]]$tfmode)
  rep_mirna <- rep(mirna_significativo, length(targets_ledge))
  
  df_union <- as.data.frame(rep_mirna)
  df_union$Target <- targets_ledge
  colnames(df_union)[1] <- "MIRNA"
  
  if (!exists("dataset_mirnas")){
    dataset_mirnas <- df_union
  } else {
    
    temp_dataset <- df_union
    
    ###MAL, HAY QUE HACER MATCH
    dataset_mirnas <- rbind(dataset_mirnas,temp_dataset)
    rm(temp_dataset)
    
  }
}

######HACER MATCH PARA OBTENER INFORMACION MUTUA

match_mrs <- merge(dataset_mirnas, network, 
                   by.x=c("MIRNA","Target"),
                   by.y=c("Regulator","Target"))

#########HACER MERGE CON EL FC DE CADA GEN
###DATOS DE MRNA EN DESEQ2 RLOG
mrna_deseq2 <-  "/Users/daniel/OneDrive\ -\ UNIVERSIDAD\ NACIONAL\ AUTÓNOMA\ DE\ MÉXICO/owncloud/X99/TCGA-BRCA/DIFF\ EXP\ MRNA\ BRCA\ NORMALES\ VS\ TUMORALES/diffexpr-results_DESEQ2.csv" 
datos_mrna_deseq2 <- read.csv(mrna_deseq2, header = TRUE)

datos_diffexp <- datos_mrna_deseq2[,c(2,4)]

mrs_foldchange <- merge(match_mrs, datos_diffexp,
                        by.x="Target",
                        by.y="Gene")

write.csv(mrs_foldchange,"DATASET_MIRNAS_SIGNIF_COOL.csv")
#


######
###cool!
#ledge_extraccion_2 <- names(ledge_extraccion_1[[mirnas_significativo]]$tfmode)


#######

mrs_signif <- filter(test, FDR < 0.05)

write.csv(mrs_signif,"MRS_SIGNIF_05.csv")

###plot_significativos

#TOP20
plot(mrs_ledge_2, as.character(mrs_signif$Regulon)
     ,cex = 2, include = c("activity"),
     pval = as.numeric(mrs_signif$FDR))

###Plot con funcion A en script externo
###FUNCIONA COOL
a(mrs_ledge_2, as.character(mrs_signif$Regulon)
     ,cex = 2, include = c("activity"),
     pval = as.numeric(mrs_signif$FDR))

######


#plot(mrs_ledge_2, 45,cex = 10)





####Quedarse solo con los mirnas con p.value < 0.05

#mirnas_signif <- filter(test, p.value < 0.05)
mirnas_signif <- as.character(mrs_signif$Regulon)

###De cada mirna extraer los genes BLANCO DE NETWORK
dataset_mirnas <- filter(network, Regulator %in% mirnas_signif)

#########HACER MERGE CON EL FC DE CADA GEN
###DATOS DE MRNA EN DESEQ2 RLOG
mrna_deseq2 <-  "/Users/daniel/OneDrive\ -\ UNIVERSIDAD\ NACIONAL\ AUTÓNOMA\ DE\ MÉXICO/owncloud/X99/TCGA-BRCA/DIFF\ EXP\ MRNA\ BRCA\ NORMALES\ VS\ TUMORALES/diffexpr-results_DESEQ2.csv" 
datos_mrna_deseq2 <- read.csv(mrna_deseq2, header = TRUE)

datos_diffexp <- datos_mrna_deseq2[,c(2,4)]

###FIX GENES SEPT

sep_genes <- dataset_mirnas

sep_genes$Target <- gsub("sep-","SEPT", as.character(sep_genes$Target),
                         fixed = TRUE)
sep_genes$Target <- gsub("SEPT0","SEPT", as.character(sep_genes$Target),
                         fixed = TRUE)
sep_genes$Target <- gsub("ago-01","AGO1", as.character(sep_genes$Target),
                         fixed = TRUE)


#test <- dataset_mirnas[grep("sep",dataset_mirnas$Target),]

mrs_foldchange <- merge(sep_genes, datos_diffexp,
                        by.x="Target",
                        by.y="Gene")


write.csv(mrs_foldchange,"DATASET_MIRNAS_SIGNIF.csv")
