#calculate the correlation between celltype ratios and clinical parameters.
#Barplot for blood cell ratio correlation analysis
library(Seurat)
library(corrplot)
library(ggplot2)

blood_all_filter <- readRDS(file= path to blood_all_cell.rds)
blood_gdt <- subset(blood_all_filter,ccg_celltype_2_main_3 == 'gdT')
blood_nk <- subset(blood_all_filter,ccg_celltype_2_main_3 == 'NK')
blood_cd4 <- subset(blood_all_filter,ccg_celltype_2_main_3 == 'CD4+ T')
blood_cd8 <- subset(blood_all_filter,ccg_celltype_2_main_3 == 'CD8+ T')
blood_mye <- subset(blood_all_filter,ccg_celltype_2_main_3 == 'Myeloid')
blood_b <- subset(blood_all_filter,ccg_celltype_2_main_3 == 'B cell')

clinical <- read.table("path to /clinical_blood.txt",head=T,sep="\t",check.names=F,stringsAsFactors = FALSE)
sc_list <- list(blood_gdt,blood_nk,blood_cd4,blood_cd8,blood_mye,blood_b)
patient <- as.character(unique(clinical[,1]))
result_c_list <- list()
result_p_list <- list()
for(m in 1:length(sc_list)){
	cell_type=as.character(unique(sc_list[[m]]@meta.data$ccg_celltype_2))
	
	celltype_frac=c()
	for(i in 1:length(patient)){
		index=which(substr(sc_list[[m]]@meta.data[,1],1,nchar(patient[i]))==patient[i])
		patient_meta=sc_list[[m]]@meta.data[index,]
		patient_celltype_frac=patient[i]
		for(j in 1:length(cell_type)){
			index_celltype=which(patient_meta$ccg_celltype_2==cell_type[j])
			length(index_celltype)/length(index)
			patient_celltype_frac=cbind(patient_celltype_frac,length(index_celltype)/length(index))
		}
		celltype_frac=rbind(celltype_frac,patient_celltype_frac)
	}
	colnames(celltype_frac)=c("patient_id",cell_type)
	celltype_frac=as.data.frame(celltype_frac)
	for(i in 2:ncol(celltype_frac)){
		celltype_frac[,i]=as.numeric(as.character(celltype_frac[,i]))
	}

	result_c_all=c()
	result_p_all=c()
	for(i in c(4,5,6,8:18)){
		quota=as.numeric(clinical[,i])
		result_c=c()
		result_p=c()
		for(j in 2:ncol(celltype_frac)){
			cor=cor.test(quota,celltype_frac[,j],method="spearman")
			c=cor$estimate
			p=cor$p.value
			result_c=cbind(result_c,c)
			result_p=c(result_p,p)
		}
		result_c_all=rbind(result_c_all,result_c)
		result_p_all=rbind(result_p_all,result_p)
	}

	colnames(result_c_all)=cell_type
	colnames(result_p_all)=cell_type
	rownames(result_c_all)=colnames(clinical)[c(4,5,6,8:18)]
	rownames(result_p_all)=colnames(clinical)[c(4,5,6,8:18)]
	result_c_list[[m]] <- result_c_all
	result_p_list[[m]] <- result_p_all
}

result_c_list <- result_c_list[!sapply(result_c_list,is.null)]
result_p_list <- result_p_list[!sapply(result_p_list,is.null)]

result_c_all <- do.call(cbind,result_c_list)
result_p_all <- do.call(cbind,result_p_list)


#calculate GZMK+ to GZMB+ CD8 T ratio correlation with clinical parameters
blood_cd8_gzm <- subset(blood_cd8,ccg_celltype_2 %in% c("GZMB+ CD8 CTL","ZNF683+ GZMB+ CD8 CTL","RGS1+ GZMK+ CD8 CTL"))
blood_cd8_gzm$celltype_main <- 'GZMB+ CD8'
blood_cd8_gzm$celltype_main[which(blood_cd8_gzm$ccg_celltype_2 %in% c("RGS1+ GZMK+ CD8 CTL"))] <- 'GZMK+ CD8'
sc_obj <- blood_cd8_gzm
sc_obj$new_cluster <- sc_obj$celltype_main

ratio <- c()
for(i in 1:length(patient)){
	index=which(substr(sc_obj@meta.data[,1],1,nchar(patient[i]))==patient[i])
	patient_meta=sc_obj@meta.data[index,'new_cluster']
	ratio_kb <- sum(patient_meta=='GZMK+ CD8')/sum(patient_meta=='GZMB+ CD8')
	ratio <- c(ratio,ratio_kb)
}

result_c=c()
result_p=c()
for(i in c(4,5,6,8:18)){
	quota=as.numeric(clinical[,i])
	cor=cor.test(quota,ratio,method="spearman")
	c=cor$estimate
	p=cor$p.value
	result_c=c(result_c,c)
	result_p=c(result_p,p)
}

result_c <- data.frame(GZMK_divide_GZMB = result_c,row.names=colnames(clinical)[c(4,5,6,8:18)],stringsAsFactors = FALSE)
result_p <- data.frame(GZMK_divide_GZMB = result_p,row.names=colnames(clinical)[c(4,5,6,8:18)],stringsAsFactors = FALSE)

result_c_all <- cbind(result_c_all,result_c)
result_p_all <- cbind(result_p_all,result_p)

result_cp_all <- data.frame(cor_r=c(t(result_c_all['Phases',])),cor_p=c(t(result_p_all['Phases',])),celltype=colnames(result_c_all),row.names=colnames(result_c_all),stringsAsFactors=FALSE)
#result_cp_all <- data.frame(cor_r=c(t(result_c_all['PGE2(pg/ml)',])),cor_p=c(t(result_p_all['PGE2(pg/ml)',])),celltype=colnames(result_c_all),row.names=colnames(result_c_all),stringsAsFactors=FALSE)
#result_cp_all <- data.frame(cor_r=c(t(result_c_all['ALT(U/L)',])),cor_p=c(t(result_p_all['ALT(U/L)',])),celltype=colnames(result_c_all),row.names=colnames(result_c_all),stringsAsFactors=FALSE)

result_cp_all$direction <- 1
result_cp_all$direction[which(result_cp_all$cor_r >= 0)] <- 'positive'
result_cp_all$direction[which(result_cp_all$cor_r < 0)] <- 'negative'

result_cp_all <- result_cp_all[order(result_cp_all$cor_r,decreasing = T), ]
result_cp_all$celltype <- factor(result_cp_all$celltype,levels=result_cp_all$celltype)
#result_cp_all$direction <- factor(result_cp_all$direction,levels=c('positive','negative'))

plot2<-ggplot(data=result_cp_all, mapping=aes(x = celltype, y = cor_r,fill = direction))+
  geom_bar(stat="identity",position=position_dodge(0.75),width=0.6)+
  coord_cartesian(ylim=c(-1,1))+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values =c("#BD81B0",'#3B95C3'))+#,"#A1D99B", "#C7EAE5", "#EB7649", "#FDD0A2", "#F4A42B",  "#42B59B", "#CC4F58"
  theme_bar_1()+
  ylab("PGE2 correlation")+
  theme(axis.text.x=element_text(angle=90, hjust=1))#+geom_text(aes(label = Percent_ratio), vjust = -0.5)
  
file_path <- 'path to result save dir /ALT_celltype_fre_cor_barplot_direction_blood.pdf'
ggsave(filename = file_path, plot = plot2, width = 26, height = 20, units = 'cm',limitsize = FALSE)


