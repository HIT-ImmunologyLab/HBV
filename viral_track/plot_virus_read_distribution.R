suppressMessages(library(ggplot2))


args <- commandArgs(T)
cov_file <- args[1]
result_file <- args[2]

data_for_plot <- read.delim(cov_file,header = F,sep = "\t")
colname_list <- c("genome_id","location","count")
colnames(data_for_plot) <- colname_list
ggplot(data = data_for_plot,aes(x=location,y=count/1000,fill="#CA4B48"))+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  geom_area(position = position_dodge(0.9),fill="#CA4B48")+
  ylab("Coverage(1000 reads)")+
  xlab("Viral segment(bp)")+
  ylim(0,1)
ggsave(result_file,width = 6,height = 3)
