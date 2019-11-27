
library(matrixStats)
library(dplyr)
library(ComplexHeatmap)
library(circlize)


datos <- read.csv("DF_FINAL_COSMIC.csv", row.names = 2)
datos <- datos[,2:length(colnames(datos))]
rownames(datos) <- gsub("hsa-","",rownames(datos),fixed = TRUE)


datos_heat <- as.matrix(datos[,1:7])

library(RColorBrewer)
library(ggplot2)

colores <- brewer.pal(9, "YlGn")

heatmap_1 <- Heatmap(datos_heat, cluster_columns = FALSE,
                     cluster_rows = TRUE,row_names_gp = gpar(fontsize = 6),
                     show_row_names = TRUE,col = colores,#colorRamp2(c(0, 15, 30), c("darkblue", "white", "firebrick1")),
                     name = "log2FoldChange" )


###Obtener order de los rows despues de hacer el clustering
r_order <- row_order(heatmap_1)
vals <- datos$index
vals_index <- datos$index[r_order]


######BARPLOT HALLMARKS
vals_hallmarks <- datos$hallmark#[r_order]


#####
heatmap_1

#row_ha = rowAnnotation(foo2 = runif(10), bar2 = anno_barplot(runif(10)))


bg <- "deepskyblue1"

###OLD CON HALLMARKS
#ha <- HeatmapAnnotation(foo = anno_empty(border = TRUE, height = unit(6, "cm")),
#                        Hallmark = anno_barplot(vals_hallmarks,gp=gpar(fill=bg),
#                                                axis_param = list(
#                                                  gp = gpar(fontsize = 5)
#                                                )),
#                        which = "row",
#                        show_annotation_name = c("Hallmark" = FALSE))

####NEW SIN HALLMARKS

ha <- HeatmapAnnotation(foo = anno_empty(border = TRUE, height = unit(6, "cm")),
                        which = "row")

top_anno <- c("","","","            Oncogenic Effect","   ","      Tumor suppressor Effect","")

top <- HeatmapAnnotation(foo2 = anno_text(top_anno, location = 0.8, just = "center",
                                          rot = 360,gp = gpar(fontsize = 8)),
                         top = colnames(datos[,1:7]), which = "column",
                         col = list(top = c("Oncogene_TSG" = "white",
                                            "Oncogene_all" = "white", 
                                            "TSG_all" = "white", 
                                            "Oncogene_Up" = "#FF000040",
                                            "TSG_Down" = "#FF000040",
                                            "Oncogene_Down" = "#8192FB",
                                            "TSG_Up" = "#8192FB")),
                         show_legend = FALSE,
                         show_annotation_name = FALSE
                         )
###########


heatmap_1 <- Heatmap(datos_heat, cluster_columns = FALSE,
                     cluster_rows = TRUE,#row_names_gp = gpar(fontsize = 6),
                     show_row_names = TRUE,top_annotation = top,col = colores,#colorRamp2(c(0, 15, 30), c("darkblue", "white", "firebrick1")),
                     name = "# Genes",
                     column_names_rot = 45 ,right_annotation = ha,
                     cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                       grid.text(datos_heat[i, j], x, y, gp = gpar(fontsize = 8))
                     },
                     row_names_gp = gpar(fontsize = 8,col = ifelse(vals> 0, "red", ifelse(vals== 0, "black", "blue"))))

ht = draw(heatmap_1)

#value = runif(10)
decorate_annotation("foo", {
  # value on x-axis is always 1:ncol(mat)
  x = 1:length(vals)#1:10
  # while values on y-axis is the value after column reordering
  value = vals_index#value[co]
  pushViewport(viewport(yscale = c(0.4, (length(vals))+0.5), xscale = c(-1.1, 1.1)))
  grid.lines(c(0, 0), c(0.4, (length(vals))+0.5), gp = gpar(lty = 3),
             default.units = "native")
  grid.points(value, rev(x), pch = 16, size = unit(2, "mm"),
              gp = gpar(col = ifelse(value > 0, "red", ifelse(value == 0, "black", "blue"))
  ), default.units = "native")
  grid.xaxis(at = c(-1, 0, 1), gp = gpar(fontsize = 5))
  grid.text("Oncogenic net activity", unit(5, "mm"),unit(-5, "mm") ,just = "center",
            gp = gpar(fontsize = 7))
  popViewport()
})

decorate_annotation("top", {
  
  #grid.text("Hallmark", unit(5, "mm"),unit(-5.5, "mm") ,just = "center",
  #          gp = gpar(fontsize = 5))
  grid.text("", unit(5, "mm"),unit(-5.5, "mm") ,just = "center",
            gp = gpar(fontsize = 5))


})



