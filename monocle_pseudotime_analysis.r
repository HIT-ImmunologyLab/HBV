library(patchwork)
library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle)
library(RColorBrewer)
library(DDRTree)
library(ggplot2)

##reading data from Seurat object
data <- as(as.matrix(seurat@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
         
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)         
# Filtering of genes detected in less than 1% of cells with a minimum expression threshold of 0.1
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(Biobase::fData(monocle_cds), num_cells_expressed > nrow(pd) * 0.05))

HSMM <- monocle_cds[expressed_genes,]


#Trajectory step 1 Screening for genes for differentiation analysis

disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)


#Trajectory step 2: reduce data dimensionality  
HSMM <- reduceDimension(HSMM, max_components = 2,method = 'DDRTree')
#Trajectory step 3: order cells along the trajectory  
HSMM <- orderCells(HSMM)
dir <- paste0("./HBV/monocle2/CD8/hbvCTL2/")
if(!file.exists(dir)){
	dir.create(dir,recursive = TRUE)
}
###############saving pictures
file <- paste0(dir,"/","group.pdf")
p <- plot_cell_trajectory(HSMM, color_by = "group_name",cell_size = 0.5)+ 
	scale_color_manual(values = c('#1BB4B8','#FF6666','#A07DB6','#79A72C','#999999'))+
	theme(legend.title = element_blank())+
	guides(color = guide_legend(override.aes = list(size = 2)))
ggsave(file,p,width =5,height = 4)


p <- plot_cell_trajectory(HSMM, show_cell_names = F, color_by = "Pseudotime") + scale_color_viridis_c()
file <- paste0(dir,"/","Pseudotime4.pdf")
ggsave(file,p,width =8,height = 6)


p <- plot_cell_trajectory(HSMM, color_by = "State",label = T)
file <- paste0(dir,"/","State.pdf")
ggsave(file,p,width =8,height = 6)

p <- plot_cell_trajectory(HSMM, color_by = "donor_name",label = T)+scale_color_manual(values=color.vector)
file <- paste0(dir,"/","donor_name.pdf")
ggsave(file,p,width =8,height = 6)

p <- plot_cell_trajectory(HSMM, color_by = "new_cluster",label = T)
file <- paste0(dir,"/","new_cluster.pdf")
ggsave(file,p,width =8,height = 6)

