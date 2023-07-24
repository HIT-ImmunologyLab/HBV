# Survival analysis for bulk RNA-seq samples
#process the microarray data
library(affy)
library(hgu133plus2.db)
library(Biobase)
library(BisqueRNA)
library(survival)
library("survminer")

dir_cels <- "path to .cel files/GSE14520_gz"
affyData <- ReadAffy(celfile.path = dir_cels)
CLLrma=mas5(affyData)
sampleNames(CLLrma) <- gsub(".CEL.gz$", "", sampleNames(CLLrma))
eset_f=exprs(CLLrma)

ids=toTable(hgu133plus2SYMBOL)
rownames(ids) <- ids[,1]
eset_f <- eset_f[intersect(rownames(eset_f),ids[,1]),]
eset_f <- as.data.frame(eset_f)
eset_f$gene_symbol <- ids[rownames(eset_f),2]
eset_f <- aggregate(x = eset_f[,-434],by = list(eset_f$gene_symbol), FUN = mean)
rownames(eset_f) <- eset_f[,1]
eset_f <- eset_f[,-1]

liver_all_final <- readRDS(file='file to liver all cell rds file')
liver_all_count <- seurat_matrix_count(liver_all_final)

eset_f <- as.matrix(eset_f)
bulk.eset <- Biobase::ExpressionSet(assayData = eset_f)
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names= colnames(liver_all_count),
                       SubjectName= liver_all_final@meta.data$donor_name,
                       cellType= liver_all_final@meta.data$celltype_final)
sc.meta <- data.frame(labelDescription=c("SubjectName","cellType"),row.names=c("SubjectName","cellType"))
sc.pdata <- new("AnnotatedDataFrame",data=sc.pheno,varMetadata=sc.meta)
sc.eset <- Biobase::ExpressionSet(assayData=liver_all_count,phenoData=sc.pdata)
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=FALSE)
estimates <- res$bulk.props

clinical <- read.table('path to bulk sample clinical information/GSE14520_Extra_Supplement.txt',header=TRUE,sep='\t',stringsAsFactors=FALSE)
rownames(clinical) <- clinical[,3]
clinical <- clinical[colnames(estimates),]

estimates_cd8 <- estimates[c('CD8_C01','CD8_C02','CD8_C03','CD8_C04','CD8_C05','CD8_C06','CD8_C07','CD8_C08','CD8_C09','CD8_C10'),]
estimates_cd4 <- estimates[c('CD4_C01','CD4_C02','CD4_C03','CD4_C04','CD4_C05','CD4_C06','CD4_C07','CD4_C08','CD4_C09','CD4_C10','CD4_C11'),]
estimates_mye <- estimates[c('Mye_C01','Mye_C02','Mye_C03','Mye_C04','Mye_C05','Mye_C06','Mye_C07','Mye_C08','Mye_C09','Mye_C10','Mye_C11','Mye_C12'),]
estimates_b <- estimates[c('Plasma B cell','Naive B cell','Memory B cell','unclear B'),]
estimates_hep <- estimates[c('Endo_1','Endo_2','Endo_3','Endo_4','Ep_1','Ep_2','Ep_3','Hep_1','Hep_2','Myofib'),]
estimates_nk <- estimates[c('cNK','IFNG+ LrNK','LrNK','PTGDS+ cNK'),]

#for example, survival analysis for CD8 T cell ratio with sample survival
estimates_subratio <- sweep(estimates_cd8,2,colSums(estimates_cd8),"/")
celltypes <- c('CD8_C01','CD8_C02','CD8_C03','CD8_C04','CD8_C05','CD8_C06','CD8_C07','CD8_C08','CD8_C09','CD8_C10')

clinical <- read.table('/home/project/HBV_result_chencg/blood/result_08_14/visualization/marker_gene/clinical/GSE14520_hcc_bulk/GSE14520_Extra_Supplement.txt',header=TRUE,sep='\t',stringsAsFactors=FALSE)
rownames(clinical) <- clinical[,3]
clinical <- clinical[colnames(estimates_subratio),]

clinical_hbv_hcc <- subset(clinical,HBV.viral.status %in% c('CC','AVR-CC'))
estimates_hbc_hcc <- estimates_subratio[,rownames(clinical_hbv_hcc)]

for(i in 1:length(celltypes)){
	#生存
	propor <- data.frame(fraction=estimates_hbc_hcc[celltypes[i],],status=clinical_hbv_hcc[colnames(estimates_hbc_hcc),'Survival.status'],row.names=colnames(estimates_hbc_hcc))
	propor$fraction[which(is.na(propor$fraction))] <- 0
	#propor$status <- factor(propor$status,levels=c(0,1))
	propor$level <- 'high'
	propor$level[which(propor$fraction < median(propor$fraction))] <- 'low'
	propor$month <- clinical_hbv_hcc[rownames(propor),'Survival.months']

	fit <- survfit(Surv(month,status) ~ level,data=propor)

	pdf(paste0("/home/project/HBV_result_chencg/blood/result_08_14/visualization/marker_gene/clinical/GSE14520_hcc_bulk/survival/hbv_hcc_status_",celltypes[i],"_Survival_median_subratio.pdf"),onefile = F)
	p <- ggsurvplot(fit,pval=TRUE)+labs(title = paste0(celltypes[i]))#+ylab('survival')
	print(p)
	dev.off()
	
	propor$level <- 'high'
	propor$level[which(propor$fraction < mean(propor$fraction))] <- 'low'
	propor$month <- clinical_hbv_hcc[rownames(propor),'Survival.months']

	fit <- survfit(Surv(month,status) ~ level,data=propor)

	pdf(paste0("path to result dir/hbv_hcc_status_",celltypes[i],"_Survival_mean_subratio.pdf"),onefile = F)
	p <- ggsurvplot(fit,pval=TRUE)+labs(title = paste0(celltypes[i]))
	print(p)
	dev.off()
}


seurat_matrix_count <- function(seurat_object){
  data_dgCMatrix <- seurat_object@assays$RNA@counts
  row_pos <- data_dgCMatrix@i+1
  col_pos <- findInterval(seq(data_dgCMatrix@x)-1,data_dgCMatrix@p[-1])+1
  val <- data_dgCMatrix@x
  
  ## conver sparse matrix to dense matrix
  data <- matrix(data=0L, nrow = data_dgCMatrix@Dim[1], ncol = data_dgCMatrix@Dim[2])
  for (i in seq_along(val)){
    data[row_pos[i],col_pos[i]] <- val[i]
  }
  row.names(data) <- data_dgCMatrix@Dimnames[[1]]
  colnames(data) <- data_dgCMatrix@Dimnames[[2]]
  return(data)
  }

  fraction=function(matrix){
  a=length(which(matrix>0))
  b=length(matrix)
  round(a/b,2)
}
