library("survival")
library("survminer")
library("dplyr")


datos <- read.csv("DATOS_PROCESADOS.csv", header = TRUE)

###Remover pacientes con tumor stage X y not_reported

datos_filtrados <- filter(datos,tumor_stage != "not reported")
datos_filtrados <- filter(datos_filtrados,tumor_stage != "stage x")
###22 pacientes sin datos de stage

datos <- datos_filtrados

datos$tumor_stage <- as.character(datos$tumor_stage)

tumor_stage <- as.character(datos$tumor_stage)

tumor_stage <- gsub("stage iv","4",tumor_stage,fixed = TRUE)
tumor_stage <- gsub("stage iiic","3",tumor_stage,fixed = TRUE)
tumor_stage <- gsub("stage iiib","3",tumor_stage,fixed = TRUE)
tumor_stage <- gsub("stage iiia","3",tumor_stage,fixed = TRUE)
tumor_stage <- gsub("stage iii","3",tumor_stage,fixed = TRUE)
tumor_stage <- gsub("stage iib","2",tumor_stage,fixed = TRUE)
tumor_stage <- gsub("stage iia","2",tumor_stage,fixed = TRUE)
tumor_stage <- gsub("stage ii","2",tumor_stage,fixed = TRUE)
tumor_stage <- gsub("stage ib","1",tumor_stage,fixed = TRUE)
tumor_stage <- gsub("stage ia","1",tumor_stage,fixed = TRUE)
tumor_stage <- gsub("stage i",1,tumor_stage,fixed = TRUE)

tumor_stage <- as.numeric(tumor_stage)

datos$stage_bien <- tumor_stage














#datos_completos <- read.csv("COMPLETOS_DATOS_CLINICOS.csv", header = TRUE)

datos$anos_bien <- ((datos$dias_bien))

colnames_bien <- gsub(".csv","", colnames(datos), fixed = TRUE)
colnames_bien <- gsub("MIR","miR", colnames_bien, fixed = TRUE)
colnames_bien <- gsub("5P","5p", colnames_bien, fixed = TRUE)
colnames_bien <- gsub("3P","3p", colnames_bien, fixed = TRUE)

colnames(datos) <- colnames_bien







#####Dataset_final

datos_completos_age <- datos


###Remover pacientes sin age at diagnosis
datos_completos_age$edad_bien <- ((datos_completos_age$age_at_diagnosis)/365)

#datos_filtrados <- datos[!is.na(datos_completos_age$edad_bien),]

#datos_completos_age <- datos_filtrados


##HACER CUT CON EL PROMEDIO DE LA EDAD

datos_completos_age$edad_bien_label <- cut(datos_completos_age$edad_bien, breaks = 2, 
           labels = c("down","up"))





#######KAPLAN MEIERS

# Survival curves


status_bien <- gsub("Alive",0,datos_completos_age$status_vital)
status_bien <- gsub("Dead",1,status_bien)

datos_completos_age$Status <- as.numeric(status_bien)



###REPETIR TODOS LOS MULTIVARIADOS, 
mirnas <- colnames(datos_completos_age)[grep("hsa",colnames(datos_completos_age))]

mirnas_bien <- mirnas[-grep("valor_expresion",mirnas)]

mirnas <- mirnas_bien

covariates <- c(mirnas,"edad_bien_label","stage_bien")


##############UNIVARIADOS
######TODOS LOS STAGES
###

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(anos_bien, Status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = datos_completos_age)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
resultado_univariado_cox_todos_stages <- as.data.frame(res)

write.csv(resultado_univariado_cox_todos_stages,"RES_COX_UNI.csv")


#res.cox <- coxph(Surv(anos_bien, Status) ~ hsa.miR.1307.3p,
#                 data = datos_completos_age)




mirnas <- colnames(datos_completos_age)[grep("hsa",colnames(datos_completos_age))]

###MULTIVARIADOS SOLO CON EDAD y stage por miRNA

#res.cox <- coxph(Surv(anos_bien, Status) ~ + edad_bien_label 
#                 + stage_bien + hsa.miR.3677.3p, 
#                 data = datos_completos_age)


covariates_multi <- covariates[grep("hsa", covariates)]
multi_formulas <- sapply(covariates_multi,
                        function(x) as.formula(paste('Surv(anos_bien, Status)~ 
                                                     edad_bien_label + 
                                                     stage_bien +', x)))
multi_models <- lapply( multi_formulas, function(x){coxph(x, data = datos_completos_age)})

multivariado_results <- lapply(multi_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res_multi <- t(as.data.frame(multivariado_results, check.names = FALSE))
resultado_multivariado_cox_todos_stages <- as.data.frame(res_multi)

######HACER UNO POR UNO

#hsa.miR.940 +
#hsa.miR.1307.3p +
#hsa.miR.340.3p +
#hsa.miR.877.5p
#falta el mir 151

res_cox_multi <- res.cox <- coxph(Surv(anos_bien, Status) ~ +edad_bien_label 
                                  + stage_bien + 
                                    hsa.miR.877.5p,
                                  data = datos_completos_age)

summary_cox <- summary(res.cox)

cox_coef <- summary_cox$coefficients

cox_intervalos <- summary_cox$conf.int[,c(3,4)]


res_cox_final <- cbind(cox_coef, cox_intervalos)
res_cox_final <- res_cox_final[,c(2,5,6,7)]
res_cox_final <- as.data.frame(res_cox_final)
res_cox_final$`lower .95` <- substr(res_cox_final$`lower .95`,1,4)
res_cox_final$`upper .95` <- substr(res_cox_final$`upper .95`,1,4)
res_cox_final$coef <- substr(res_cox_final$`exp(coef)`,1,4)
res_cox_final$HR <- paste(res_cox_final$`lower .95`,res_cox_final$`upper .95`, sep="-")
res_cox_final$HR <- paste0("(",res_cox_final$HR,")")


res_cox_final$HR <- paste(res_cox_final$coef,res_cox_final$HR, sep=" ")

res_cox_final <- res_cox_final[,c(6,2)]
colnames(res_cox_final)[2] <- "p.value"
res_cox_final$p.value <- substr(res_cox_final$p.value,1,5)

res_cox_final
summary(res.cox)


#############
###########
res.cox <- coxph(Surv(anos_bien, Status) ~ +edad_bien_label 
                 + stage_bien + 
                   hsa.miR.151a.5p +
                   hsa.miR.940 +
                   hsa.miR.1307.3p +
                   hsa.miR.340.3p +
                   hsa.miR.877.5p
                 , 
                 data = datos_completos_age)

summary_cox <- summary(res.cox)

cox_coef <- summary_cox$coefficients

cox_intervalos <- summary_cox$conf.int[,c(3,4)]


res_cox_final <- cbind(cox_coef, cox_intervalos)
res_cox_final <- res_cox_final[,c(2,5,6,7)]
res_cox_final <- as.data.frame(res_cox_final)
res_cox_final$`lower .95` <- substr(res_cox_final$`lower .95`,1,4)
res_cox_final$`upper .95` <- substr(res_cox_final$`upper .95`,1,4)
res_cox_final$coef <- substr(res_cox_final$`exp(coef)`,1,4)
res_cox_final$HR <- paste(res_cox_final$`lower .95`,res_cox_final$`upper .95`, sep="-")
res_cox_final$HR <- paste0("(",res_cox_final$HR,")")


res_cox_final$HR <- paste(res_cox_final$coef,res_cox_final$HR, sep=" ")

res_cox_final <- res_cox_final[,c(6,2)]
colnames(res_cox_final)[2] <- "p.value"
res_cox_final$p.value <- substr(res_cox_final$p.value,1,5)


write.csv(res_cox_final,"RES_COX_MULTI.csv")


