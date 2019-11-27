
library(dplyr)


stage1 <- read.csv("normales_vs_stage1.csv", header = TRUE)

stage2 <- read.csv("normales_vs_stage2.csv", header = TRUE)

stage3 <- read.csv("normales_vs_stage3.csv", header = TRUE)

stage4 <- read.csv("normales_vs_stage4.csv", header = TRUE)


#######Filtros, mirna con foldchange de al menos 1 y padj de 0.05

filtro_stage1 <- filter(stage1, log2FoldChange > 0.5 & padj < 0.05)

filtro_stage2 <- filter(stage2, log2FoldChange > 0.5 & padj < 0.05)

filtro_stage3 <- filter(stage3, log2FoldChange > 0.5 & padj < 0.05)

filtro_stage4 <- filter(stage4, log2FoldChange > 0.5 & padj < 0.05)


########ENCONTRAR COMPARTIDOS


compartidos <- Reduce(function(x, y) merge(x, y, by.x = "X", by.y = "X"), list(filtro_stage1, 
                                                          filtro_stage2, 
                                                          filtro_stage3,
                                                          filtro_stage4))



write.csv(compartidos, "RESULTADOS_MIRNAS_SOBRE_FC_0_5.csv")




