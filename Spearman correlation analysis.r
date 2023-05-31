###############
library(patchwork)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(Seurat)
library(ggrepel)

#######read Read blood scTCR data
clone_meta <- read.table("Liver_PBMC_combine_rename_id.txt",sep="\t",header=T,
                        stringsAsFactor=F)
clone_meta$tcr_full_length <- paste0(clone_meta$chain1_V,clone_meta$chain1_J,clone_meta$chain1_C,clone_meta$chain1_CDR3,
                                  clone_meta$chain2_V,clone_meta$chain2_J,clone_meta$chain2_C,clone_meta$chain2_CDR3)

clone_sub <- clone_meta%>%filter(loc=="PBMC")%>%filter(majorCluster%in%CD8_celltype)
CD8_celltype <- c("GZMB+ CD8 CTL","ZNF683+ GZMB+ CD8 CTL","RGS1+ GZMK+ CD8 CTL","naive CD8","CD8 Tem")


#######Read the experimentally acquired HBV-specific  TCR data
specific_tcr <- read.table("./Tetramer/TCR_result.txt",sep="\t",header=T)
specific_tcr$tcr_full_length <- paste0(specific_tcr$alpha_v_gene,specific_tcr$alpha_j_gene,specific_tcr$alpha_c_gene,specific_tcr$alpha_cdr3,
specific_tcr$beta_v_gene,specific_tcr$beta_j_gene,specific_tcr$beta_c_gene,specific_tcr$beta_cdr3)


clone_sub$HBV_specific <- "non_specific"
clone_sub$HBV_specific[clone_sub$tcr_full_length%in%specific_tcr$tcr_full_length]="HBV_specific"
blood_specific <- clone_sub%>%filter(HBV_specific=="HBV_specific")


#######plot pictures
p <- list(NULL)
for(sample in c(unique(blood_specific$sample))){
	
	blood_specific_sample <- blood_specific%>%filter(patient==sample)
	clonetype_index <- data.frame(tcr_full_length=table(blood_specific_sample$tcr_full_length)%>%names(),
									count=table(blood_specific_sample$tcr_full_length)%>%as.character()%>%as.numeric())
	clonetype_index$rank <- rank(-c(clonetype_index$count),ties.method="min")
	clonetype_index$group <- "scrna_data"
	specific_tcr_sub <- specific_tcr%>%filter(tcr_full_length%in%blood_specific_sample$tcr_full_length)
	expriment_result <- data.frame(tcr_full_length=table(specific_tcr_sub$tcr_full_length)%>%names(),
									count=table(specific_tcr_sub$tcr_full_length)%>%as.character()%>%as.numeric())
	expriment_result$rank <- rank(-c(expriment_result$count),ties.method="min")
	expriment_result$group <- "expriment_data"
	dat <- left_join(clonetype_index,expriment_result,by="tcr_full_length")
	xlab_name <- "count(scrna-seq)"
	ylab_name <- "count(expriment)"
	p[[sample]]<- ggscatter(dat, x = "count.x", y = "count.y", 
			  repel = TRUE,#Avoid label overlap
			  color = "black", size = 0.8, # Set the color and size of the dot
			  add = "reg.line",  # Add regression line
			  add.params = list(color = "black", fill = "lightgray"), # The color of the regression line is set to red, and the color of the interval is set to gray
			  conf.int = TRUE, # Add confidence interval for regression line
			  cor.coef = TRUE, # Add correlation coefficient
			  cor.coeff.args = list(method = "spearman",label.sep = "\n"),cor.coef.size = 3 ,cor.coef.color = "red")+
			 labs(x=xlab_name,y=ylab_name)+
			  ggtitle(paste0("donor: ",sample))+
			  theme(plot.title = element_text(face = "bold",family = "Times"))+
			  scale_y_continuous(breaks=c(0,1,5,10,20,50,100,500,1000), trans="log1p")+
			  scale_x_continuous(breaks=c(0,1,5,10,20,50,100,500), trans="log1p")

}

p <- p[!sapply(p,is.null)]
merge_plot=grid.arrange(grobs = p, ncol = 3, nrow = 3)
ggsave("count.pdf",merge_plot,width=12,height=9)


		