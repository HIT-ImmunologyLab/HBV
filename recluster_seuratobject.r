library(Seurat)
library(harmony)


recluster_fun=function(seurat_object,dim_list,res_list,result_dir,prefix){
  recluster_object=seurat_object
  dim_list=dim_list
  res_list=res_list
  result_dir=result_dir
  prefix=prefix
  
  recluster_object=NormalizeData(recluster_object, assay=NULL, 
  normalization.method = "LogNormalize", scale.factor = 10000)
  
  ## select variable gene
  seurat_object_subset <- FindVariableFeatures(recluster_object,selection.method = "vst",nfeatures=2000)
  ## scale data

  seurat_object_subset <- ScaleData(seurat_object_subset,vars.to.regress = "percent.mt",verbose = TRUE)
  ## PCA
  seurat_object_subset <- RunPCA(seurat_object_subset,pc.genes=VariableFeatures(object = seurat_object_subset),verbose = TRUE)
  for(i in 1:length(dim_list)){
  dim_value=dim_list[i]
  dim_seurat_object_subset=RunHarmony(seurat_object_subset,"orig.ident",plot_convergence = TRUE)
  dim_seurat_object_subset <- RunUMAP(dim_seurat_object_subset, dims = 1:dim_value,reduction = "harmony")    
  dim_seurat_object_subset <-FindNeighbors(dim_seurat_object_subset, dims = 1:dim_value,reduction = "harmony")
  

  for(j in 1:length(res_list)){
  resolution_value=res_list[j]
  dim_seurat_object_subset <-FindClusters(dim_seurat_object_subset,resolution = resolution_value)
  identity(dim_seurat_object_subset)
  
save_dir=paste(result_dir,"/","dim_value",dim_value,"_res",resolution_value,sep="")
dir.create(save_dir)
saveRDS(dim_seurat_object_subset,paste(save_dir,"/",prefix,"_umap_seurat_object.rds",sep=""))






png(paste(save_dir,"/","umap_cluster.png",sep=""),height = 4000,width=6000,res=300)
a=UMAPPlot(object = dim_seurat_object_subset, pt.size = 0.51, label = TRUE)    #umap 
print(a)
dev.off()


png(paste(save_dir,"/","umap_donor.png",sep=""),height = 4000,width=6000,res=300)
a=UMAPPlot(object = dim_seurat_object_subset, pt.size = 0.51, label = TRUE,group.by="orig.ident")    #group by samples
print(a)
dev.off()


png(paste(save_dir,"/","umap_group.png",sep=""),height = 4000,width=6000,res=300)
a=UMAPPlot(object = dim_seurat_object_subset, pt.size = 0.51, label = TRUE,group.by="group_name")    #group by phase
print(a)
dev.off()



png(paste(save_dir,"/","umap_split_donor.png",sep=""),height = 1000,width=30000,res=150)
a=UMAPPlot(object = dim_seurat_object_subset, pt.size = 0.51, label = TRUE,split.by="orig.ident")    #split by samples
print(a)
dev.off()



png(paste(save_dir,"/","umap_split_group.png",sep=""),height = 2000,width=12000,res=150)
a=UMAPPlot(object = dim_seurat_object_subset, pt.size = 0.51, label = TRUE,split.by="group_name")    #split by phase
print(a)
dev.off()

att_gene_list=unique(c("CD3E","CD3D","CD8A","CD4","LYZ","CD3D","CD3E","CD3G","CD247","NCAM1","FCGR3A","NFIL3","CD19",
"GZMA","GZMB","GZMK","GZMH","GNLY","PRF1","FCGR3A","S1PR1","S1PR2","S1PR3","S1PR4","S1PR5","KLF2","PRDM1","ITGAE","ITGA1",
"PDCD1","CTLA4","LAG3","TOX","TCF7","TIGIT","HAVCR2","XCL1","RGS1","EOMES","TBX21","CCR7","LEF1","SELL","RORC","ICAM1",
"VCAM1","CCL20","ICAM1","ZNF683","CX3CR1","CXCR6","ITGAX","FCGR3B","PTPRC",
"CD3D","CD3E","CD8A","NKG7","NCAM1","SLC4A10","TRDV2",
"CD79A","MS4A1","CD19","CST3","CD68","LYZ","VCAN","FCN1","ST14","CLEC10A",
"CD1C","BATF3","FCER1A","CD68","MAFB","CSF1R","CD14","ITGAM","FCGR3A","GP9","PF4","ALB","CD27","SDC1","IGHG1",
"IGHA2","FOXP3","CD4","CTLA4","CXCR6","CX3CR1","IL7R","TNFRSF4","CCR7","SELL","TCF7","LEF1","KLRC1","IFNG"))
att_gene_list=intersect(att_gene_list,rownames(dim_seurat_object_subset))


png(paste(save_dir,"/","FeaturePlot.png",sep=""),height = 30000,width=30000,res=300)
a=FeaturePlot(dim_seurat_object_subset, features = att_gene_list,reduction="umap")
print(a)
dev.off()



seurat_object_subset_markers <- FindAllMarkers(dim_seurat_object_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(seurat_object_subset_markers,paste(save_dir,"/","cluster_markers.txt",sep=""),
sep="\t",row.names=F,quote=F)
  
  }
  
	}
}
##################


result_dir="***/seurat_object/Liver"###result object save directory
dim_list=c(50,40,30)##Select the dimensions you want to output
res_list=c(0.5,1)##Select the resolutions you want to output
prefix="all_liver"##Record the source of your sample
recluster_fun(combine_object,dim_list,res_list,result_dir,prefix)




