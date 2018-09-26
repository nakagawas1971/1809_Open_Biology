##############################################
##0. create sh script for importing bed files#
##############################################
library (data.table)
setwd("~/GoogleDrive/Share/review_paper/R_analyses")

#import metadata
meta <- fread("metadata.tsv")

#download bed file from ENCODE server
meta_peaks <- meta[meta$`File format` == "bed narrowPeak"
              & meta$Assembly == "GRCh38",]
download <- paste("wget",meta_peaks$`File download URL`)
write.table(download1,file="download_bed.txt",sep="\t",row.names=F,col.names=F,quote=F) #download using Terminal sh command




#############################################################
##1. create data table containing NEAT1 bed fold change data#
#############################################################
library (data.table)
setwd("~/GoogleDrive/Share/review_paper/R_analyses/bed")
fl <- list.files(pattern = ".bed")
NEAT1_start <- 65422798
NEAT1_end <- 65445538
NEAT1_bed <- data.table(position = NEAT1_start:NEAT1_end) #create data table containing genome position

#convert bed file format into data table containing eCLIP FC value for each base
for (i in 1:length(fl)) {
  print(paste(i,"/",length(fl)))
  ID <- gsub("\\..+$", "", fl[i]) #extract ID name from file name
  
  #read bed file
  temp <- fread(fl[i]) 
  
  #select peak value corresponding to NEAT1 region
  temp_NEAT1 <- temp[temp$V1=="chr11"&temp$V2>NEAT1_start&temp$V3<NEAT1_end,]
  
  #prepare column for each RBP
  eval(parse(text = paste("NEAT1_bed$",ID," <- 0",sep = "")))
  
  #put bed fc value into data table
  for (j in 1:nrow(temp_NEAT1)){
    start <- grep(temp_NEAT1$V2[j],NEAT1_bed$position)
    end <- grep(temp_NEAT1$V3[j],NEAT1_bed$position)
    NEAT1_bed[start:end,(i+1)] <- temp_NEAT1$V7[j]
  }
}
write.table(NEAT1_bed,file="../NEAT1_bedcount.txt",sep="\t",row.names=F,col.names=T,quote=F)



##############################################################################
#2. create data table containing average eCLIP FC value avarage for each RBPs#
#   write coverage plot for each RBPs                                        #
##############################################################################
library(ggplot2)
library(data.table)
setwd("~/GoogleDrive/Share/review_paper/R_analyses")

#import metadata
meta <- fread("metadata.tsv")

#get list of RBP name
meta_HepG2 <- meta_peaks[meta_peaks$`Biosample term name` == "HepG2"]
RBP_HepG2 <- sub("-human","",unique(meta_HepG2$`Experiment target`)) # list of RBP for HepG2 eCLIP

meta_K562 <- meta_peaks[meta_peaks$`Biosample term name` == "K562"]
RBP_K562 <- sub("-human","",unique(meta_K562$`Experiment target`)) # list of RBP for K562 eCLIP

#HepG2
#plot % of input and prepare average value table
NEAT1_start <- 65422798
NEAT1_end <- 65445538
NEAT1_bed <- fread("NEAT1_bedcount.txt")
ncol(NEAT1_bed)
HepG2_bed <- data.table(position = NEAT1_start:NEAT1_end) #prepare data table for fc
rep1_rn <- grep(1, meta_HepG2$`Biological replicate(s)`) #collect rn for eCLIP rep1 
rep2_rn <- grep(2, meta_HepG2$`Biological replicate(s)`) #collect rn for eCLIP rep2

for (i in 1:length(RBP_HepG2)){
  print(paste("HepG2",i,"/",length(RBP_HepG2))) #print log
  RBP_temp <- RBP_HepG2[i] #pick up the name of RBPs
  IDs <- meta_HepG2$`File accession`[grep(RBP_temp,meta_HepG2$`Experiment target`)] # pick up the experiment IDs for each RBP_temp
  RBP_temp_rn <- grep(RBP_temp,sub("-human","",meta_HepG2$`Experiment target`)) # pike up row numbers for each RBP_temp
  
  if (RBP_temp == "HNRNPU"){ # exceptional case (grep command cannnot discriminate HNRNPU and HNRNPUL1)
    RBP_temp_rn <- c(125, 126) #HNRNPU rn = 125 and 126
  }
  
  File_ID_1 <- meta_HepG2$`File accession`[intersect(RBP_temp_rn,rep1_rn)] #obtain file accession ID for eCLUP rep1
  File_ID_2 <- meta_HepG2$`File accession`[intersect(RBP_temp_rn,rep2_rn)] #obtain file accession ID for eCLUP rep2

  eval(parse(text = paste("value1 <- NEAT1_bed$",File_ID_1,sep=""))) #extract data for eCLIP rep 1
  eval(parse(text = paste("value2 <- NEAT1_bed$",File_ID_2,sep=""))) #extract data for eCLIP rep 2
  

    value <- (value1+value2)/2 #get average
    value[value < 0] <- 0 #adjust minus FC value to 0
    eval(parse(text = paste("HepG2_bed$",RBP_temp," <- value",sep = ""))) # create table containing average eCLIP FC value
    
    g <- ggplot(NULL)
    g <- g + geom_line (data = HepG2_bed,aes(x = HepG2_bed$position, y = value), size=0.4, color = "coral")
    
    g = g +
      theme(axis.text.y= element_text(size=8),axis.title.y= element_text(size=9)) +
      theme(axis.text.x= element_text(size=8),axis.title.x= element_text(size=9)) +
      ylab("log2(eCLIP/input)") +
      xlab("") +
      labs(title = paste(RBP_HepG2[i]," (HepG2)",sep = "")) +
      theme(plot.title = element_text(size=9)) +
      ylim(0,7)
    g
    ggsave(filename = paste("./png/NEAT1/HepG2_bed_",RBP_HepG2[i],".png",sep = ""),width = 8, height = 2)
}

write.table(HepG2_bed,file="NEAT1_HepG2_bed.txt",sep="\t",row.names=F,col.names=T,quote=F)

#K562
#plot % of input and prepare average value table
NEAT1_start <- 65422798
NEAT1_end <- 65445538
NEAT1_bed <- fread("NEAT1_bedcount.txt")
ncol(NEAT1_bed)
K562_bed <- data.table(position = NEAT1_start:NEAT1_end) #prepare data table for fc
rep1_rn <- grep(1, meta_K562$`Biological replicate(s)`) #collect rn for eCLIP rep1 
rep2_rn <- grep(2, meta_K562$`Biological replicate(s)`) #collect rn for eCLIP rep2

for (i in 1:length(RBP_K562)){
  print(paste("K562",i,"/",length(RBP_K562))) #print log
  RBP_temp <- RBP_K562[i] #pick up the name of RBPs
  IDs <- meta_K562$`File accession`[grep(RBP_temp,meta_K562$`Experiment target`)] # pick up the experiment IDs for each RBP_temp
  RBP_temp_rn <- grep(RBP_temp,sub("-human","",meta_K562$`Experiment target`)) # pike up row numbers for each RBP_temp
  
  if (RBP_temp == "HNRNPU"){ # exceptional case (grep command cannnot discriminate HNRNPU and HNRNPUL1)
    RBP_temp_rn <- c(125, 126) #HNRNPU rn = 125 and 126
  }
  
  File_ID_1 <- meta_K562$`File accession`[intersect(RBP_temp_rn,rep1_rn)] #obtain file accession ID for eCLUP rep1
  File_ID_2 <- meta_K562$`File accession`[intersect(RBP_temp_rn,rep2_rn)] #obtain file accession ID for eCLUP rep2
  
  eval(parse(text = paste("value1 <- NEAT1_bed$",File_ID_1,sep=""))) #extract data for eCLIP rep 1
  eval(parse(text = paste("value2 <- NEAT1_bed$",File_ID_2,sep=""))) #extract data for eCLIP rep 2
  
  
  value <- (value1+value2)/2 #get average
  value[value < 0] <- 0 #adjust minus FC value to 0
  eval(parse(text = paste("K562_bed$",RBP_temp," <- value",sep = ""))) # create table containing average eCLIP FC value
  
  g <- ggplot(NULL)
  g <- g + geom_line (data = K562_bed,aes(x = K562_bed$position, y = value), size=0.4, color = "coral")
  
  g = g +
    theme(axis.text.y= element_text(size=8),axis.title.y= element_text(size=9)) +
    theme(axis.text.x= element_text(size=8),axis.title.x= element_text(size=9)) +
    ylab("log2(eCLIP/input)") +
    xlab("") +
    labs(title = paste(RBP_K562[i]," (K562)",sep = "")) +
    theme(plot.title = element_text(size=9)) +
    ylim(0,7)
  g
  ggsave(filename = paste("./png/NEAT1/K562_bed_",RBP_K562[i],".png",sep = ""),width = 8, height = 2)
}

write.table(K562_bed,file="NEAT1_K562_bed.txt",sep="\t",row.names=F,col.names=T,quote=F)
