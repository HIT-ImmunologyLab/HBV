##########计算clone diversity曲线
suppressMessages(library(ggrepel))
suppressMessages(library(alakazam))
suppressMessages(library(shazam))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(igraph))
suppressMessages(library(scoper))
suppressMessages(library(RColorBrewer))
suppressMessages(library(reshape2))
library(ggsci)
plotDiversityCurve_correct=function (data, colors = NULL, main_title = "Diversity", legend_title = "Group", 
log_x = FALSE, log_y = FALSE, xlim = NULL, ylim = NULL, annotate = c("none", 
	"depth"), score = c("diversity", "evenness"), silent = FALSE, 
...) 
{
annotate <- match.arg(annotate)
score <- match.arg(score)
if (all(is.na(data@groups)) || length(data@groups) == 1) {
	group_labels <- NA
}
else if (annotate == "none") {
	group_labels <- setNames(data@groups, data@groups)
}
else if (annotate == "depth") {
	group_labels <- setNames(paste0(data@groups, " (n=", 
		data@n, ")"), data@groups)
}
if (score == "diversity") {
	y_value <- "d"
	y_min <- "d_lower"
	y_max <- "d_upper"
	y_label <- expression(""^q * D)
}
else if (score == "evenness") {
	y_value <- "e"
	y_min <- "e_lower"
	y_max <- "e_upper"
	y_label <- expression(""^q * e)
}
.x <- NULL
if (!all(is.na(group_labels))) {
	p1 <- ggplot(data@diversity, aes_string(x = "q", y = y_value, 
		group = data@group_by)) + ggtitle(main_title) + baseTheme() + 
		xlab("q") + ylab(y_label)  + 
		geom_line(aes_string(color = data@group_by))
	if (!is.null(colors)) {
		p1 <- p1 + scale_color_manual(name = legend_title, 
			labels = group_labels, values = colors) + scale_fill_manual(name = legend_title, 
			labels = group_labels, values = colors)
	}
	else {
		p1 <- p1 + scale_color_discrete(name = legend_title, 
			labels = group_labels) + scale_fill_discrete(name = legend_title, 
			labels = group_labels)
	}
}
else {
	if (!is.null(colors) & length(colors) == 1) {
		line_color <- colors
	}
	else {
		line_color <- "black"
	}
	p1 <- ggplot(data@diversity, aes_string(x = "q", y = y_value)) + 
		ggtitle(main_title) + baseTheme() + xlab("q") + ylab(y_label) + 
		geom_line(color = line_color)
}
if (log_x) {
	p1 <- p1 + scale_x_continuous(trans = scales::log2_trans(), 
		limits = xlim, breaks = scales::trans_breaks("log2", 
			function(x) 2^x), labels = scales::trans_format("log2", 
			scales::math_format(2^.x)))
}
else {
	p1 <- p1 + scale_x_continuous(limits = xlim)
}
if (log_y) {
	p1 <- p1 + scale_y_continuous(trans = scales::log2_trans(), 
		limits = ylim, breaks = scales::trans_breaks("log2", 
			function(x) 2^x), labels = scales::trans_format("log2", 
			scales::math_format(2^.x)))
}
else {
	p1 <- p1 + scale_y_continuous(limits = ylim)
}
p1 <- p1 + do.call(theme, list(...))
if (!silent) {
	plot(p1)
}
invisible(p1)
}



###################liver hbvctl tcr diversity
immune_info_data = read.table("liver_ctl_meta.txt",sep="\t",header=T,stringsAsFactor=F)

result_dir <- "~/clone_abundance_analysis"
result_prefix <- "Liver_TCR"
AHB_colors <- rep("#E41A1C",3)
IA_colors <- rep("#377EB8",3)
IC_colors <- rep("#4DAF4A",4)
IT_colors <- rep("#FF7F00",4)

patient_colors <- c(AHB_colors,IA_colors,IC_colors,IT_colors)

## clonal abundance
immune_info_data[,'patient'] <- immune_info_data[,'donor_name']
immune_info_data[,'patient'] <- apply(as.matrix(immune_info_data[,'patient']), 1, function(x){
unlist(strsplit(x,"_"))[1]
})



## clonal diversity
diversity_curve_immune_patient <- alphaDiversity(immune_info_data, group="patient",
							  min_q=0, max_q=20, step_q=0.1,ci=0.95, nboot=200)


diversity_curve_immune_patient@diversity$patient=factor(diversity_curve_immune_patient@diversity$patient,
levels=c("AHB01" ,"AHB02","AHB03","A04","A06","A07","C01","C03","C04","C06","T02", "T04","T05","T06"))
p1= plotDiversityCurve_correct(diversity_curve_immune_patient, colors = patient_colors, legend_title="patient")+guides(colour =guide_legend(override.aes=list(size=5))) +
theme(plot.title=element_text(hjust=0.5,size=12))+ggtitle("HBV-specific TCR clonal diversity(liver)")




###################blood hbvctl tcr diversity


immune_info_data=read.table("blood_ctl_meta.txt",sep="\t",header=T,stringsAsFactor=F)
result_dir <- "~/clone_abundance_analysis"
result_prefix <- "Blood_TCR"

AHB_colors <- rep("#E41A1C",3)
IA_colors <- rep("#377EB8",3)
IC_colors <- rep("#4DAF4A",5)
IT_colors <- rep("#FF7F00",4)

patient_colors <- c(AHB_colors,IA_colors,IC_colors,IT_colors)



## clonal abundance
immune_info_data[,'donor_name'] <- immune_info_data[,'patient']
immune_info_data[,'patient'] <- apply(as.matrix(immune_info_data[,'patient']), 1, function(x){
unlist(strsplit(x,"_"))[1]
})


## clonal diversity
diversity_curve_immune_patient2 <- alphaDiversity(immune_info_data, group="patient",
								  min_q=0, max_q=20, step_q=0.1,ci=0.95, nboot=200)


diversity_curve_immune_patient2@diversity$patient=factor(diversity_curve_immune_patient2@diversity$patient,
levels=c("AHB01" ,"AHB02","AHB03","A04","A06","A07","C01","C03","C04","C05","C06","T02", "T04","T05","T06"))




save_file <- paste(result_dir,"/",result_prefix,"_clonal_diversity_by_patientxx.pdf",sep="")
p2= plotDiversityCurve_correct(diversity_curve_immune_patient2, colors = patient_colors, legend_title="patient")+
guides(colour =guide_legend(override.aes=list(size=5))) +
theme(plot.title=element_text(hjust=0.5,size=12))+ggtitle("HBV-specific TCR clonal diversity(blood)")

p=p1/p2
ggsave(save_file,p,width=5,height=10)
