library(LymphoSeq)
library(dplyr)

LymphoSeq::clonalRelatedness
#################################blood
clone_meta = read.table("E:\\desktop\\tcr\\blood_ctl_meta.txt",sep="\t",header=T,
                        stringsAsFactor=F)
result = list()
for(patient in unique(clone_meta$patient_id)){
  print(patient)
  a = clone_meta[clone_meta$patient_id%in%patient,]
  
  a1=a[!grepl("\\|C",a$chain2_CDR3),]
  a2=a[grep("\\|C",a$chain2_CDR3),]
  if(nrow(a2)>0){
    temp_a = list()
    for(i in 1:nrow(a2)){
      strings=strsplit(a2[i,"chain2_CDR3"],"\\|")%>%unlist()
      temp1=a2[i,]
      temp1$chain2_CDR3 = strings[1]
      temp2=a2[i,]
      temp2$chain2_CDR3 = strings[2]
      temp_a[[i]]=rbind(temp1,temp2)
    }
    temp_a=do.call(rbind,temp_a)
    
    a=rbind(a1,temp_a)
  }
  b = group_by(a,chain2_CDR3,,sampleType)%>%summarise(count=n())
  b$nucleotide=b$chain2_CDR3
  
  result[[patient]] = b
}


clonalRelatedness=function (list, editDistance = 3) 
{
  
  relatedness <- as.numeric()
  clonality <- as.numeric()
  cloneResolved <- as.character()
  i <- 1
  for (i in 1:length(list)) {
    file <- list[[i]]
    top <- file[order(file$count, decreasing = TRUE), "nucleotide"]%>%as.data.frame()%>%.[1,1]
    result1 <- length(which(stringdist::stringdist(top, file$nucleotide, 
                                                   method = "lv") < editDistance))/nrow(file)
    relatedness <- c(relatedness, result1)
  }
  df <- data.frame(samples = names(list), clonalRelatedness = relatedness)
  return(df)
  
}
aa = clonalRelatedness(list = result, editDistance = 3)




colnames(clone_meta)[16]="samples"
relatedness_value = left_join(aa,clone_meta[,c("samples","sampleType")],by="samples")
relatedness_value_blood = relatedness_value[!duplicated(relatedness_value),]


clone_meta = read.table("E:\\desktop\\tcr\\ctl_meta.txt",sep="\t",header=T,
                        stringsAsFactor=F)
result = list()
for(patient in unique(clone_meta$patient_id)){
  print(patient)
  a = clone_meta[clone_meta$patient_id%in%patient,]
  
  a1=a[!grepl("\\|C",a$chain2_CDR3),]
  a2=a[grep("\\|C",a$chain2_CDR3),]
  if(nrow(a2)>0){
    temp_a = list()
    for(i in 1:nrow(a2)){
      strings=strsplit(a2[i,"chain2_CDR3"],"\\|")%>%unlist()
      temp1=a2[i,]
      temp1$chain2_CDR3 = strings[1]
      temp2=a2[i,]
      temp2$chain2_CDR3 = strings[2]
      temp_a[[i]]=rbind(temp1,temp2)
    }
    temp_a=do.call(rbind,temp_a)
    
    a=rbind(a1,temp_a)
  }
  b = group_by(a,chain2_CDR3,,sampleType)%>%summarise(count=n())
  b$nucleotide=b$chain2_CDR3
  
  result[[patient]] = b
}


clonalRelatedness=function (list, editDistance = 3) 
{
  
  relatedness <- as.numeric()
  clonality <- as.numeric()
  cloneResolved <- as.character()
  i <- 1
  for (i in 1:length(list)) {
    file <- list[[i]]
    top <- file[order(file$count, decreasing = TRUE), "nucleotide"]%>%as.data.frame()%>%.[1,1]
    result1 <- length(which(stringdist::stringdist(top, file$nucleotide, 
                                                   method = "lv") < editDistance))/nrow(file)
    relatedness <- c(relatedness, result1)
  }
  df <- data.frame(samples = names(list), clonalRelatedness = relatedness)
  return(df)
  
}
aa = clonalRelatedness(list = result, editDistance = 3)




colnames(clone_meta)[16]="samples"
relatedness_value = left_join(aa,clone_meta[,c("samples","sampleType")],by="samples")
relatedness_value_blood = relatedness_value[!duplicated(relatedness_value),]



#####################################
####################################
###################################liver

clone_meta = read.table("E:\\desktop\\tcr\\ctl_meta.txt",sep="\t",header=T,
                        stringsAsFactor=F)
result = list()
for(patient in unique(clone_meta$donor_name)){
  print(patient)
  a = clone_meta[clone_meta$donor_name%in%patient,]
  
  a1=a[!grepl("\\|C",a$chain2_CDR3),]
  a2=a[grep("\\|C",a$chain2_CDR3),]
  if(nrow(a2)>0){
    temp_a = list()
    for(i in 1:nrow(a2)){
      strings=strsplit(a2[i,"chain2_CDR3"],"\\|")%>%unlist()
      temp1=a2[i,]
      temp1$chain2_CDR3 = strings[1]
      temp2=a2[i,]
      temp2$chain2_CDR3 = strings[2]
      temp_a[[i]]=rbind(temp1,temp2)
    }
    temp_a=do.call(rbind,temp_a)
    
    a=rbind(a1,temp_a)
  }
  b = group_by(a,chain2_CDR3,,sampleType)%>%summarise(count=n())
  b$nucleotide=b$chain2_CDR3
  
  result[[patient]] = b
}


clonalRelatedness=function (list, editDistance = 3) 
{
  
  relatedness <- as.numeric()
  clonality <- as.numeric()
  cloneResolved <- as.character()
  i <- 1
  for (i in 1:length(list)) {
    file <- list[[i]]
    top <- file[order(file$count, decreasing = TRUE), "nucleotide"]%>%as.data.frame()%>%.[1,1]
    result1 <- length(which(stringdist::stringdist(top, file$nucleotide, 
                                                   method = "lv") < editDistance))/nrow(file)
    relatedness <- c(relatedness, result1)
  }
  df <- data.frame(samples = names(list), clonalRelatedness = relatedness)
  return(df)
  
}
aa = clonalRelatedness(list = result, editDistance = 3)




colnames(clone_meta)[1]="samples"
relatedness_value = left_join(aa,clone_meta[,c("samples","sampleType")],by="samples")
relatedness_value_liver = relatedness_value[!duplicated(relatedness_value),]%>%filter(samples!="A02_LIVER_ZLY")%>%filter(samples!="C05_LIVER_MFB")

relatedness_value_liver$pool="Liver hbvctl TCR pool"
relatedness_value_liver$samples=sapply(relatedness_value_liver$samples,function(x){strsplit(x,split = "_")%>%unlist()%>%.[1]})

relatedness_value_blood$pool="Blood hbvctl TCR pool"


relatedness_data=rbind(relatedness_value_liver,relatedness_value_blood)


                                                      
library(ggsci)
library(ggrepel)


p = ggplot(relatedness_data,aes(x=sampleType, y=clonalRelatedness,color=sampleType))+geom_boxplot(size=0.5)+
  geom_jitter(aes(x=sampleType, y=clonalRelatedness,fill=sampleType),size =1.5,shape = 21,, width = 0.3,stroke = 0.1)+
  #geom_text_repel(aes(x=sampleType,y=clonalRelatedness,label=samples),color="black",size=1)+
  theme_classic()+facet_wrap(.~pool,scales = "free_x",nrow = 1)+
  theme(strip.text=element_text(size=12,face="bold"),
        strip.background = element_rect(fill = 'white', colour = 'white', size = rel(2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=4),
        legend.key.height=unit(0.5,"cm"),
        axis.text=element_text(size=7),
  )+
  scale_fill_uchicago(palette = "dark")+
  xlab("")+
  scale_color_uchicago(palette = "dark")+scale_y_log10()
ggsave("E:\\desktop\\pool_ctl.pdf",p,width=7,height=3)



