#source("http://bioconductor.org/biocLite.R")
#biocLite("Heatplus")
#install.packages('gplots')

library(gplots)
library(Heatplus)
library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)
library (data.table)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(ggdendro)

setwd("~/GoogleDrive/Share/review_paper/R_analyses")


#read eCLIP bed fold change value table
NEAT1_HepG2_bed <- fread ("NEAT1_HepG2_bed.txt", header = TRUE)
NEAT1_HepG2_bed <- NEAT1_HepG2_bed[,-1] #remove column containing position information
colnames(NEAT1_HepG2_bed) <- paste(colnames(NEAT1_HepG2_bed),"_HepG2",sep = "") #obtain column names == RBP names

NEAT1_K562_bed <- fread ("NEAT1_K562_bed.txt", header = TRUE)
NEAT1_K562_bed <- NEAT1_K562_bed[,-1] #remove column containing position information
colnames(NEAT1_K562_bed) <- paste(colnames(NEAT1_K562_bed),"_K562",sep = "") #obtain column names == RBP names

PSP_list <- fread("PSPs_list.txt",header = FALSE) # import list of paraspeckle proteins

#create data table containing eCLIP FC values
NEAT1_all_bed <- cbind(NEAT1_K562_bed,NEAT1_HepG2_bed)
NEAT1_all_bed <- data.table(NEAT1_all_bed)

NEAT1_start <- 65422798
NEAT1_end <- 65445538

##################################################
#create binned sheet (bin = 50 window/25 overlap)#
##################################################

#import data
col <- colnames(NEAT1_all_bed) #obtain column names == RBP names
num_bin <- seq(1,nrow(NEAT1_all_bed),by = 25) #create bins from NEAT1 start to NEAT1_end by 25
NEAT1_bed_bin <- as.data.frame(NEAT1_all_bed[1:(length(num_bin)-1),]) #prepare data table for eCLIP FC value for each RBPs

for (i in 1:(length(num_bin)-2)){
  NEAT1_temp <- NEAT1_all_bed[num_bin[i]:num_bin[i+2],] #prepare temporal table corresponding to each bin (i+2 = 50 nt window, 25 bin)
  binsum <- apply(NEAT1_temp,2,mean) #calculate mean value (e.g. average eCLIP FC value for each bin)
  NEAT1_bed_bin[i,] <- binsum #put values into data table
}
#write data
write.table(NEAT1_bed_bin,file="NEAT1_bed_bin_50_25.txt",sep="\t",row.names=F,col.names=T,quote=F)

#prepare data table for "bound" RBPs (RBP with eCLIP FC value >log2(5))
NEAT1_bed_bin <- read.table("NEAT1_bed_bin_50_25.txt",header = TRUE)
NEAT1_bed_bin <- as.matrix(NEAT1_bed_bin)
NEAT1_bound_bed_bin <- NEAT1_bed_bin[,(apply(NEAT1_bed_bin,2,max) > log2(5))] #select Neat1-binding RBPs cutoff: log2(5)
NEAT1_bound_bed_bin <- as.matrix(NEAT1_bound_bed_bin)

#prepare transformed matrix
NEAT1_bound_bed_bin_tr <- t(NEAT1_bound_bed_bin)


#######################################################
#clustering analyses (RBP) and drow heatmap/dendrogram#
#######################################################

#RBP clustering analyses
NEAT1_bound_bed_bin_tr_dist <- dist(NEAT1_bound_bed_bin_tr) #calculate distance value
NEAT1_bound_bed_bin_tr_clust <- hclust (d = NEAT1_bound_bed_bin_tr_dist, method = "ward.D2") #clustering using dist value by ward.D2 or average
RBP_class <- cutree(tree=NEAT1_bound_bed_bin_tr_clust,k=6) #cut RBP into 6 branches
RBP_order <- rev(NEAT1_bound_bed_bin_tr_clust$order) #get RBP order number
NEAT1_bound_bed_bin_tr_sorted <- NEAT1_bound_bed_bin_tr[RBP_order,] #sort RBP with clustering number

#plot RBP dendrogram using ggdendro
ggdendrogram(NEAT1_bound_bed_bin_tr_clust, rotate = TRUE, size = 2)
ddata <- dendro_data(NEAT1_bound_bed_bin_tr_clust, type = "rectangle")
p <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(breaks=NULL) +
  scale_x_reverse(breaks=NULL) +
  #geom_tile() +
  xlab("") +
  ylab("") +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )
p
ggsave("NEAT1_RBP_dendrogram.png",width = 2, height = 18)

############################
#plot heatmap using ggplot2#
############################

#prepare vector "RBP_color" containing color list to highlight PSPs in the RBP list
RBP_list <- unique(row.names(NEAT1_bound_bed_bin_tr_sorted)) #RBP list: RBP names in heat map
RBP_list <- sub("_K562", "", RBP_list)
RBP_list <- sub("_HepG2", "", RBP_list)
RBP_color <- c(1:length(RBP_list)) #RBP_color: color vector for RBP in heatmap
RBP_color[] <- "black"
for (i in 1:nrow(PSP_list)){ #search for PSPs
  PSPs_num <- grep(PSP_list[i],RBP_list) #pick up RBP_list number matching the PSP
  if (length(PSPs_num) != 0) {
    RBP_color[PSPs_num] <- "blue"
  }
}

#prepare melted data table for ggplot2
NEAT1_bound_bed_bin_tr_sorted_melted <- melt(t(NEAT1_bound_bed_bin_tr_sorted)) #prepare melted table

#draw heatmap (tiles) using ggplot2
g <- ggplot(NEAT1_bound_bed_bin_tr_sorted_melted, aes(Var1*25, Var2) ) +
  geom_tile(aes(fill = value)) +
  theme(axis.text.x = element_text(size = 10),axis.title.x=element_text(size=12)) + # X label size
  theme(axis.text.y=element_text(size=10, color = RBP_color),axis.title.y=element_text(size=12)) + #y label size
  theme(panel.background = element_blank()) +
  scale_fill_gradientn("value", colours = rev(brewer.pal(9, "Spectral")), na.value = "white") +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("") +
  ylab("RBPs")
g
ggsave(filename = "NEAT1_RBP_heatmap.png", g, width = 10, height = 15)



##############################################################
#clustering analyses (NEAT1 bins) and drow heatmap/dendrogram#
##############################################################

#clustering analyses of NEAT1 bins 
NEAT1_bound_bed_bin_dist <- dist(NEAT1_bound_bed_bin) #calculate distance value
NEAT1_bound_bed_bin_clust <- hclust (d = NEAT1_bound_bed_bin_dist, method = "ward.D2") #clustering using dist value by ward.D2 or average
Neat1_class <- cutree(tree=NEAT1_bound_bed_bin_clust,k=4) #cut NEAR1 region into 4 branches


#draw NEAT1 dendrogram using ggdendro
ggdendrogram(NEAT1_bound_bed_bin_clust, rotate = TRUE, size = 2)
ddata <- dendro_data(NEAT1_bound_bed_bin_clust, type = "rectangle")
p <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(breaks=NULL) +
  scale_x_reverse(breaks=NULL) +
  #geom_tile() +
  xlab("") +
  ylab("") +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )
p
ggsave("NEAT1_domain_dendrogram.png",width = 4, height = 6)


#draw Neat1 eCLIP domain cut with 4 branches
g1 <-ggplot(domain, aes(x = bin*25, y = 1, fill = class)) +
  scale_fill_manual(values = c("#00CFB4", "springgreen2", "gray90", "violet"), guide = "none") +
  geom_tile() +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks=NULL) +
  scale_x_continuous(breaks=NULL) +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.grid = element_blank())

ggsave("NEAT1_eCLIP_domains.png", g1, width = 10, height = 1)


#draw repeat in NEAT1 region
NEAT1_repeat <- fread("NEAT1_repeat.txt") #repeatmasker bed file downloaded from table browser
NEAT1_start <- 65422798
NEAT1_end <- 65445538

g2 <-ggplot() +
  geom_rect(data = NEAT1_repeat, aes(xmin = NEAT1_repeat$V2 - NEAT1_start, xmax = NEAT1_repeat$V3 - NEAT1_start, ymin = 0, ymax = 1, fill = "gray")) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks=NULL) +
  scale_x_continuous(breaks=NULL) +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.grid = element_blank()) +
  theme(legend.position = 'none') 
  
g2 <- g2 +
  geom_line(data = NEAT1,aes(x= position - NEAT1_start,y= 1),color="gray40",size=1)
g2
  ggsave("NEAT1_repeats.png", g2, width = 10, height = 0.5)

