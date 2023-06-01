
#Validation of zonated gene signatures using scRNA-seq (edit by ccg)
library(Seurat)
library(scales)
squash_axis <- function(from, to, factor) { 
    # Args:
    #   from: left end of the axis
    #   to: right end of the axis
    #   factor: the compression factor of the range [from, to]

  trans <- function(x) {    
      # get indices for the relevant regions
      isq <- x > from & x < to
      ito <- x >= to

      # apply transformation
      x[isq] <- from + (x[isq] - from)/factor
      x[ito] <- from + (to - from)/factor + (x[ito] - to)

      return(x)
  }

  inv <- function(x) {
      # get indices for the relevant regions
      isq <- x > from & x < from + (to - from)/factor
      ito <- x >= from + (to - from)/factor

      # apply transformation
      x[isq] <- from + (x[isq] - from) * factor
      x[ito] <- to + (x[ito] - (from + (to - from)/factor))

      return(x)
  }

# return the transformation
  return(trans_new("squash_axis", trans, inv))
}

seurat <- readRDS(file='path to liver rds/liver_final.rds')
seurat[["central_liver"]] <- PercentageFeatureSet(seurat, features = c('GLUL','CYP2E1'))
seurat[["portal_liver"]] <- PercentageFeatureSet(seurat, features = c('ASS1','ASL','ALB'))
seurat[["Xenobiotic metabolism"]] <- PercentageFeatureSet(seurat, features = c('ACSM2B','BAAT','CES2','CYP1A2','CYP2E1','CYP3A4','CYP3A5','EPHX1','GSTA2','SULT2A1'))
seurat[["Fatty acid biosynthesis"]] <- PercentageFeatureSet(seurat, features = c('ACSM2B','CD74','CYP1A2','CYP2E1','CYP3A4'))
seurat[["Secretion"]] <- PercentageFeatureSet(seurat, features = c('ALB','FGA','FGB','FGG','HP','IGF1','IL1RN','ORM1','PLA2G2A','TF','TIMP1','TTR'))
seurat[["Iron homeostasis"]] <- PercentageFeatureSet(seurat, features = c('HAMP','HPX','TF'))
seurat_sub <- subset(seurat,subset=seurat_clusters %in% c(1,2,11,14))

save_file <- paste0('path to save file/liver_hyp_mark_function_featureplot.png')
pathways <- c('central_liver','portal_liver','Xenobiotic.metabolism','Fatty.acid.biosynthesis','Secretion','Iron.homeostasis','IL6ST','SAA4')
p = list(NULL)
for(i in 1:length(pathways)){
	p[[i]] <- FeaturePlot(object = seurat_sub, features = pathways[i], reduction = 'umap',label = TRUE,ncol=1,pt.size=0.1)+coord_trans(y = squash_axis(-3, 10, 10),x=squash_axis(2, 8, 10))+theme(axis.text.y = element_text(size=8,color="black"),axis.text.x = element_text(size=8,color="black"))
}
p <- p[!sapply(p,is.null)]
p <- p[[1]] + p[[2]] + p[[3]] + p[[4]] + p[[5]] + p[[6]] + p[[7]] + p[[8]] + plot_layout(ncol = 2)
ggsave(filename = save_file, plot = p, width =10, height = 12, dpi=300,limitsize = FALSE)



