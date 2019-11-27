
##ESTOS DATOS YA ESTAN NORMALIZADOS Y EXPRESADOS EN LOG2

datos <- read.csv("vst_transformation.csv",row.names = 1,
                  stringsAsFactors = FALSE)


mirnas_cool <- read.csv("RESULTADOS_MIRNAS_SOBRE_FC_0_5.csv", header = TRUE,
                        row.names = 1)

mirnas_chidos <- c(as.character(mirnas_cool$X),"grupo")

datos_filtrados <- datos[mirnas_chidos,]

rownames(datos_filtrados) <- gsub("-",".",rownames(datos_filtrados),fixed = TRUE)

datos <- as.data.frame(t(datos_filtrados),
                       stringsAsFactors = FALSE)

grupo <- datos$grupo
grupo <- gsub("normales","1",grupo)
grupo <- gsub("tumores1","2",grupo)
grupo <- gsub("tumores2","3",grupo)
grupo <- gsub("tumores3","4",grupo)
grupo <- gsub("tumores4","5",grupo)

datos$grupo <- grupo

datos_numericos <- as.data.frame(sapply(datos, as.numeric))
rownames(datos_numericos) <- rownames(datos)

#shapiro
#test_norm <- rnorm(1000, 3, .25)

shapiro.test(datos_numericos$hsa.let.7g.3p)

library(dplyr)
library(Hmisc)
#regresion <- lm(hsa.let.7a.1 ~ grupo, data = datos)

#xyplot sirve para hacer el plot
#xyplot(hsa.mir.183 ~ grupo, data = datos)
boxplot(hsa.miR.7.5p ~ grupo, data = datos_numericos,
        xlab="Stage",ylab="log2 raw read counts",
        main="hsa-miR-7-5p")
#
#cor.test(~hsa.let.7a.1 + grupo, data=datos)

###esto ya esta bien
datos_matrix <- as.matrix(datos_numericos)

matriz_correlacion <- rcorr(datos_matrix, type = "pearson")#c("pearson","spearman"))

#Extraer datos de P y R
matriz_correlacion_R <- matriz_correlacion$r
matriz_correlacion_P <- matriz_correlacion$P

#################

#solo hay que extraer la columna "grupo" de R
matriz_correlacion_R <- as.data.frame(subset(matriz_correlacion_R, select = c(grupo)))
#cambiar nombre de grupo a grupo_r
colnames(matriz_correlacion_R) <- c("grupo_r")

#solo hay que extraer la columna "grupo" de P
matriz_correlacion_P <- as.data.frame(subset(matriz_correlacion_P, select = c(grupo)))
#cambiar nombre de grupo a grupo_r
colnames(matriz_correlacion_P) <- c("grupo_Pval")

####unir dataframe de matriz correlacion R y P

matriz_final_correlacion <- cbind(matriz_correlacion_R,matriz_correlacion_P)

write.csv(matriz_final_correlacion,"matriz_final_correlacion_pearson.csv")

