### R
#### commands for IBD block analysis
# in this script, there will be several examples of plots aimed to show the IBD sharing patterns between populations.
# information for each population will be used.
# First, create a file that for each ID associates a population, some grouping categories (continent, macro greographic region, language family, etc), and geographic coordinate to plot on a map.

library(data.table)

# Read and concatenate files
folder_path<-"/Users/vborda/Documents/Public_data/INS_data/IBD_Analyses/IBD_segments/hapibd/results/"
file_list <- list.files(folder_path, pattern = "\\.ibd$", full.names = TRUE)
df_combined <- rbindlist(lapply(file_list, fread))

colnames(df_combined)<-c("firstID","firstHapIndex","secondID","secondHapIndex", "chromosome", "start","end","length_no_merged","length")
ibd<-df_combined[,-8]

path<-"/Users/vborda/Documents/Public_data/INS_data/data+RefPops/"
key<-read.table(paste0(path,"FULL_1KGP_HGDP_PGP_IDs.txt"), head=FALSE, fill = TRUE,sep = "\t")
colnames(key)<-c("ID","Population","Acronym","City")
### PEL includes related individuals.
pel_unrelated<-read.table("/Users/vborda/Documents/Public_data/1000Genomes/IDs/Lista_ID_POP_1000Genomes_PEL.txt", header = FALSE)
colnames(pel_unrelated)<-c("id","Population")

library(dplyr)

key <- key %>%
  mutate(Population = ifelse(
    Population == "PEL" & !(ID %in% pel_unrelated$id),
    "PEL_related",
    Population
  ))


regions<-c("Iquitos","Tumbes","Lambayeque","Trujillo","Ancash","PEL","Lima","Ica","Arequipa","Ayacucho","Moquegua","Tacna","Cusco","Puno")
demo.peru<-key[key$Population %in% regions,c(1,2)]
library(dplyr)
demo.peru <- demo.peru %>% mutate(Population = gsub("Ancash", "Huaraz", Population))
demo.peru <- demo.peru %>% mutate(Population = gsub("Ica", "El Carmen", Population))
demo.peru <- demo.peru %>% mutate(Population = gsub("Puno", "Juliaca", Population))

colnames(demo.peru)[1]<-"firstID"
ibdmatch<-merge(ibd,demo.peru,by.x = "firstID",by.y = "firstID")   # associate the population source for the first sample ID of the couple
colnames(ibdmatch)[9]<-"source1"
colnames(demo.peru)[1]<-"secondID"
ibdmatch2<-merge(ibdmatch,demo.peru,by.x = "secondID",by.y = "secondID") # associate the population source for the second sample ID of the couple
colnames(ibdmatch2)[10]<-"source2"

FAM_unrelated<-read.table("/Users/vborda/Documents/Public_data/INS_data/Freeze_data/2023Dec_INS_LDGH_1KGP_HGDP_autosomes_updatedIDs.fam")
FAM_unrelated<-FAM_unrelated[,c(1,2)]
colnames(FAM_unrelated)<-c("Population","id")
FAM_unrelated <- FAM_unrelated %>%
  mutate(Population = ifelse(
    Population == "PEL" & !(id %in% pel_unrelated$id),
    "PEL_related",
    Population
  ))

FAM_unrelated <- FAM_unrelated %>% mutate(Population = gsub("Ancash", "Huaraz", Population))
FAM_unrelated <- FAM_unrelated %>% mutate(Population = gsub("Afro_des", "El Carmen", Population))
FAM_unrelated <- FAM_unrelated %>% mutate(Population = gsub("Puno", "Juliaca", Population))
poporder<-c("Iquitos","Tumbes","Lambayeque","Trujillo","Huaraz","PEL","Lima","El Carmen","Arequipa","Ayacucho","Moquegua","Tacna","Cusco","Juliaca")
FAM_unrelated<-FAM_unrelated[FAM_unrelated$Population %in% poporder,]
library(dplyr)

count_data <- FAM_unrelated %>%
  group_by(Population) %>%
  summarise(count = n())

tmp<-ibdmatch2[ibdmatch2$firstID %in% FAM_unrelated$id,]
ibdmatch2<-tmp[tmp$secondID %in% FAM_unrelated$id,]
ibd_10gen<-ibdmatch2[ibdmatch2$length > 15,]              ### nine generations Generationtime = 3/2*Lenght
ibd_9gen<-ibdmatch2[ibdmatch2$length > 16.6,]            ### nine generations Generationtime = 3/2*Lenght
ibd_8gen<-ibdmatch2[ibdmatch2$length > 18.75,]          ### eigth generations Generationtime = 3/2*Lenght 
ibd_7gen<-ibdmatch2[ibdmatch2$length > 21.42,]
ibd_6gen<-ibdmatch2[ibdmatch2$length > 25,]            ### six generations Generationtime = 3/2*Lenght


ibd_older<-ibdmatch2[ibdmatch2$length < 18.75,] 
ibd_older<-ibd_older[ibd_older$length > 3,] 

### table sharing per population
ibd<-ibd_older
#ibd<-ibd_10gen
#ibd<-ibd_9gen
ibd<-ibd_8gen
#ibd<-ibd_7gen
#ibd<-ibd_6gen
poporder<-c("Iquitos","Tumbes","Lambayeque","Trujillo","Huaraz","PEL","Lima","El Carmen","Arequipa","Ayacucho","Moquegua","Tacna","Cusco","Juliaca")
#pops<-table(demo.peru$Population)
pops<-table(FAM_unrelated$Population)

perpop<-matrix(NA,length(poporder),11)

colnames(perpop)<-c("population","samplesize","numberSharingTot","numbersharingWithin","numberSharingOut","FreqSharingTot","FreqsharingWithin","FreqSharingOut","Mean_lengthsharingWithin","totallenghtsharing","howmanypops")
perpop[,1]<-poporder
perpop[,2]<-pops[poporder]

for (i in 1:nrow(perpop)){
popp<-poporder[i]
within<-which(ibd$source1%in%popp & ibd$source2%in%popp)
tempWithin<-ibd[within,]
tempTOT<-ibd[union(which(ibd$source1%in%popp),which(ibd$source2%in%popp)),]
tempOut<-tempTOT[-which(tempTOT$source1 == tempTOT$source2),]
perpop[i,3]<-nrow(tempTOT)
perpop[i,4]<-nrow(tempWithin)
perpop[i,5]<-nrow(tempOut)
perpop[i,6]<-nrow(tempTOT)/as.numeric(perpop[i,2])
perpop[i,7]<-nrow(tempWithin)/as.numeric(perpop[i,2])
perpop[i,8]<-nrow(tempOut)/as.numeric(perpop[i,2])
perpop[i,9]<-mean(tempTOT$length)
popvarie<-c(tempTOT$source1,tempTOT$source2)
perpop[i,10]<-sum(tempTOT$length)
perpop[i,11]<-length(unique(popvarie))
}

outpath<-"/Users/vborda/Documents/Public_data/INS_data/IBD_Analyses/IBD_segments/hapibd/plots_and_tables/"

write.table(perpop,paste0(outpath,"202508_hapibd_16cM_popInfoIBDsharing.txt"),sep="\t", row.names = F, quote=F)
#perpop<-read.table("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/IBD_analyses/Rephased_phasedibd_popInfoIBDsharing.txt",sep="\t",header=T, as.is=T,comment.char = "", fill=T, quote="")

#------------------------------------------------------------------
 ### SECTION 1 visualize exchange between populations as number of events
 #------------------------------------------------------------------
library(ggplot2)
 # create a file to plot a symmetric matrix of exchange
 bestmirror<-ibd
 bestmirror$source1<-ibd$source2
 bestmirror$source2<-ibd$source1
 bestdouble<-rbind(ibd,bestmirror)
 
 # now plot as symmetric matrix 
 p<-ggplot(bestdouble, aes(x = source1, y = source2,size=length)) + 
   labs(x="source1", y="source2", title="IBD sharing") +
   geom_point(shape=21) +
   scale_color_gradient(low="lightblue", high="red") +  # the color code corresponds to the LOD score
   theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                      axis.text.y=element_text(size=9),
                      plot.title=element_text(size=11),
                      legend.text=element_text(size=7)) +
   theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
   scale_x_discrete(limits=poporder)+       #the poporder can be a different population order that you want to display in the plot
   scale_y_discrete(limits=poporder)
 p
 ggsave(paste0(outpath,"2025_hapibd_7G_IBDsharingblocksSymmetricPlay.pdf"), useDingbats=FALSE)
 
# the plot above did not have the diagonal (within pop sharing)
 # this plot shows the sharing of segments within each population, color code for LOD, length on y axis
# ibd<-ibd[which(ibd$LOD>10),]
 
diagonal<-ibd[which(ibd$source1==ibd$source2),]
pintern<-ggplot(diagonal, aes(x = source1,y=length)) + #,color=LOD)) + 
  labs(x="source1", title="IBD sharing within population ") +
  geom_point(shape=21) +
  scale_color_gradient(low="lightblue", high="red") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11),
                     legend.text=element_text(size=7)) +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  scale_x_discrete(limits=poporder)
pintern
ggsave(paste0(outpath,"2025_hapibd_9G_IBDsharingblocksWithinPopPlay.pdf"), useDingbats=FALSE)

#------------------------------------------------------------------
# Matrices of exchange between populations
#------------------------------------------------------------------

# matrix with the total number of shared blocks
matrixIBD<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-ibd[union(which(ibd$source1==pop.i),which(ibd$source2==pop.i)),]
    if (pop.i==pop.k){
      matrixIBD[i,k]<- length(which(temp$source1==temp$source2))
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBD[i,k]<-nrow(tempp)
    }
  }
}
write.table(matrixIBD,paste0(outpath,"2025_hapibd_16cM_matrix_refinedIBD_merge_sharing.txt"), sep="\t")

# make a matrix with the average length of blocks
matrixIBDAverageLength<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-ibd[union(which(ibd$source1==pop.i),which(ibd$source2==pop.i)),]
      if (pop.i==pop.k){
        temp2<-temp[(which(temp$source1==temp$source2)),]
      matrixIBDAverageLength[i,k]<- mean(temp2$length)
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBDAverageLength[i,k]<-mean(tempp$length)
    }
  }
}
write.table(matrixIBDAverageLength,paste0(outpath,"2025_hapibd_16cM_matrix_IBD_averageLength.txt"), sep="\t")


# make a matrix with the TOTAL length of blocks
matrixIBDTotLength<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-ibd[union(which(ibd$source1==pop.i),which(ibd$source2==pop.i)),]
    if (pop.i==pop.k){
      temp2<-temp[(which(temp$source1==temp$source2)),]
      matrixIBDTotLength[i,k]<- sum(temp2$length)
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBDTotLength[i,k]<-sum(tempp$length)
    }
  }
}
write.table(matrixIBDTotLength,paste0(outpath,"2025_hapibd_16cM_matrix_IBD_totalLength.txt"), sep="\t")

#pops<-table(infoID$PopName)
#pops<-pops[which(pops>0)]
pops<-pops[rownames(matrixIBD)]

#adjust for population size
matrixIBDadjustpopsize<-matrixIBD
for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
      matrixIBDadjustpopsize[i,k]<- matrixIBD[i,k]/(pops[i]*pops[k])
  }
}
write.table(matrixIBDadjustpopsize,paste0(outpath,"2025_hapibd_16cM_matrix_IBDsharingAdjustPopSize.txt"), sep="\t")

matrixIBDadjustpopsizelength<-matrixIBDTotLength
for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    matrixIBDadjustpopsizelength[i,k]<- matrixIBDTotLength[i,k]/(pops[i]*pops[k])
  }
}
write.table(matrixIBDadjustpopsizelength,paste0(outpath,"2025_hapibd_16cM_matrix_IBDsharingAdjustPopSizelength.txt"), sep="\t")

# ---------------
# Melt everything for ggplot

library(reshape2)
library(ggplot2)
meltIBD<-reshape2::melt(matrixIBD)
colnames(meltIBD)<-c("source1", "source2", "n_sharing")

meltIBDaverage<-reshape2::melt(matrixIBDAverageLength)
meltIBD$averageLength<-meltIBDaverage$value
meltIBDlength<-reshape2::melt(matrixIBDTotLength)
meltIBD$totalLength<-meltIBDlength$value
meltIBDadjuxt<-reshape2::melt(matrixIBDadjustpopsize)
meltIBD$sharingadjust<-meltIBDadjuxt$value
meltIBDlengthadjuxt<-reshape2::melt(matrixIBDadjustpopsizelength)
meltIBD$lengthadjust<-meltIBDlengthadjuxt$value


meltIBD2<-meltIBD[-(which(meltIBD$source1==meltIBD$source2)),] #exclude same pop sharing

# you can filter for more significant pairs, like pairs that share more than once, or more than the median
meltIBD22<-meltIBD2[which(meltIBD2$n_sharing!=0),]
meltIBD22<-meltIBD2[which(meltIBD2$n_sharing>1),] # more than once
meltIBD22_median<-meltIBD22[which(meltIBD22$sharingadjust>median(meltIBD22$sharingadjust)),] #more than the median

gg_8G <- ggplot(meltIBD22, aes(x = source1, y = source2, fill = sharingadjust, size = lengthadjust)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low = "lightblue", high = "red") +
  theme_bw() +
  scale_x_discrete(limits = rev(poporder)) +
  scale_y_discrete(limits = rev(poporder)) +
  xlab("Population 1") +
  ylab("Population 2") +
  labs(
    fill = "Number of IBD Shared\nadjusted by\nPopulation size",
    size = "Length (cM) adjusted by\nPopulation size"
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1, size = 7.5),  # Axis tick labels
    axis.text.y  = element_text(size = 7.5),                                    # Y tick labels
    axis.title.x = element_text(size = 8.5),                                    # X-axis title
    axis.title.y = element_text(size = 8.5),                                    # Y-axis title
    legend.title = element_text(size = 7.5),                                    # Legend title
    legend.text  = element_text(size = 6.5),                                    # Legend labels
    legend.key.size = unit(0.25, "cm")                                         # Legend box size
  )



gg_older <- ggplot(meltIBD22, aes(x = source1, y = source2, fill = sharingadjust, size = lengthadjust)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low = "lightblue", high = "red") +
  theme_bw() +
  scale_x_discrete(limits = rev(poporder)) +
  scale_y_discrete(limits = rev(poporder)) +
  xlab("Population 1") +
  ylab("Population 2") +
  labs(
    fill = "Number of IBD Shared\nadjusted by\nPopulation size",
    size = "Length (cM) adjusted by\nPopulation size"
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1, size = 7.5),  # Axis tick labels
    axis.text.y  = element_text(size = 7.5),                                    # Y tick labels
    axis.title.x = element_text(size = 8.5),                                    # X-axis title
    axis.title.y = element_text(size = 8.5),                                    # Y-axis title
    legend.title = element_text(size = 7.5),                                    # Legend title
    legend.text  = element_text(size = 6.5),                                    # Legend labels
    legend.key.size = unit(0.25, "cm")                                         # Legend box size
  )



gg_older
gg_7G
gg_8G
gg_9G
gg_10G

#meltIBD22_median_6G<-meltIBD22_median
#meltIBD22_median_7G<-meltIBD22_median
meltIBD22_median_8G<-meltIBD22_median
#meltIBD22_median_9G<-meltIBD22_median
#meltIBD22_median_10G<-meltIBD22_median
meltIBD22_median_older<-meltIBD22_median

gg_median_7G<-ggplot(meltIBD22_median_7G,aes(x=source1, y=source2, fill=sharingadjust, size=lengthadjust))+
  geom_point(shape=21)+
  scale_fill_gradient(low="lightblue", high="red") +
  scale_size_continuous(guide = guide_legend(ncol = 2)) + 
  theme_bw() +
  theme(legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.33, "cm") ,
        legend.spacing = unit(0.05, "cm"),
        axis.text.x = element_text(size = 6,angle = 45, vjust=1, hjust=1),#,colour=as.character(info$color))) +
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6)) +
  ggtitle("~7 Generations ago") +
  scale_x_discrete(limits=rev(poporder))+
  scale_y_discrete(limits=rev(poporder))+ xlab("Population 1") + ylab("Population 2") + 
  labs(fill = "Number of IBD Shared\nadjusted by\nPopulation size") + labs(size = "Length (cM) nadjusted by\nPopulation size")


gg_median_8G<-ggplot(meltIBD22_median_8G,aes(x=source1, y=source2, fill=sharingadjust, size=lengthadjust))+
  geom_point(shape=21)+
  scale_fill_gradient(low="lightblue", high="red") +
  scale_size_continuous(guide = guide_legend(ncol = 2)) + 
  theme_bw() +
  theme(legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.4, "cm") ,
        legend.spacing = unit(0.05, "cm"),
        axis.text.x = element_text(size = 6,angle = 45, vjust=1, hjust=1),#,colour=as.character(info$color))) +
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6)) +
  ggtitle("After 8 Generations ago") +
  scale_x_discrete(limits=rev(poporder))+
  scale_y_discrete(limits=rev(poporder))+ xlab("Population 1") + ylab("Population 2") + 
  labs(fill = "Number of IBD Shared\nadjusted by\nPopulation size") + labs(size = "Length (cM) nadjusted by\nPopulation size")

save(gg_median_8G,file="/Users/vborda/Documents/Public_data/INS_data/Manuscript/R2_Figures/SubFigures/Figure_3_above_ggMedian_8G.RData")


gg_median_9G<-ggplot(meltIBD22_median_9G,aes(x=source1, y=source2, fill=sharingadjust, size=lengthadjust))+
  geom_point(shape=21)+
  scale_fill_gradient(low="lightblue", high="red") +
  scale_size_continuous(guide = guide_legend(ncol = 2)) + 
  theme_bw() +
  theme(legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.33, "cm") ,
        legend.spacing = unit(0.05, "cm"),
        axis.text.x = element_text(size = 6,angle = 45, vjust=1, hjust=1),#,colour=as.character(info$color))) +
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6)) +
  ggtitle("~9 Generations ago") +
  scale_x_discrete(limits=rev(poporder))+
  scale_y_discrete(limits=rev(poporder))+ xlab("Population 1") + ylab("Population 2") + 
  labs(fill = "Number of IBD Shared\nadjusted by\nPopulation size") + labs(size = "Length (cM) nadjusted by\nPopulation size")

gg_median_10G<-ggplot(meltIBD22_median_10G,aes(x=source1, y=source2, fill=sharingadjust, size=lengthadjust))+
  geom_point(shape=21)+
  scale_fill_gradient(low="lightblue", high="red") +
  scale_size_continuous(guide = guide_legend(ncol = 2)) +  # First legend in 2 columns
  #scale_fill_gradientn(colours = terrain.colors(7))+
  # scale_fill_distiller(palette = "RdPu")+
  theme_bw() +
  theme(legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.33, "cm") ,
        legend.spacing = unit(0.05, "cm"),
        axis.text.x = element_text(size = 6,angle = 45, vjust=1, hjust=1),#,colour=as.character(info$color))) +
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6)) +
  ggtitle("After 10 Generations") +
  scale_x_discrete(limits=rev(poporder))+
  scale_y_discrete(limits=rev(poporder))+ xlab("Population 1") + ylab("Population 2") + 
  labs(fill = "Number of IBD Shared\nadjusted by\nPopulation size") + labs(size = "Length (cM) adjusted by\nPopulation size")

gg_median_older<-ggplot(meltIBD22_median_older,aes(x=source1, y=source2, fill=sharingadjust, size=lengthadjust))+
  geom_point(shape=21)+
  scale_fill_gradient(low="lightblue", high="red") +
  scale_size_continuous(guide = guide_legend(ncol = 2)) +  # First legend in 2 columns
  #scale_fill_gradientn(colours = terrain.colors(7))+
  # scale_fill_distiller(palette = "RdPu")+
  theme_bw() +
  theme(legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.4, "cm") ,
        legend.spacing = unit(0.05, "cm"),
        axis.text.x = element_text(size = 6,angle = 45, vjust=1, hjust=1),#,colour=as.character(info$color))) +
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6)) +
  ggtitle("Before 8 Generations ago") +
  scale_x_discrete(limits=rev(poporder))+
  scale_y_discrete(limits=rev(poporder))+ xlab("Population 1") + ylab("Population 2") + 
  labs(fill = "Number of IBD Shared\nadjusted by\nPopulation size") + labs(size = "Length (cM) adjusted by\nPopulation size")

save(gg_median_older,file="/Users/vborda/Documents/Public_data/INS_data/Manuscript/R2_Figures/SubFigures/Figure_3_above_ggMedian_older.RData")


gg_median_older

gg_median_10G
gg_median_9G
gg_median_8G
gg_median_7G
gg_median_6G

library(patchwork)




library(cowplot)

# Combine the plots, ensuring the legend is shared
combined_plot <- plot_grid(
  gg_median_7G, gg_median_8G, gg_median_9G, gg_median_10G,
  ncol = 2,  # Set number of columns
  rel_widths = c(1, 1),  # Adjust relative widths if necessary
  labels = c("A", "B", "C", "D")  # Optional: add labels to your plots
)# Print the combined plot
combined_plot





library(ggplot2)
library(hrbrthemes)
library(tidyverse)
library(Hmisc)
library(reshape2)
library(ggpubr)

#ggarrange(gg,ggplot() + theme_void(),gg_median,                      # First row with scatter plot
#          labels = c("A","", "B"),                  # Second row with box and dot plots
#          ncol = 1, heights = c(1, 0.1, 1)       # Labels of the scatter plot
#)

#output_path<-"/Users/vborda/Documents/Public_data/INS_data/Manuscript/Figures"
#ggsave(paste0(output_path,"/Fig_S4_IBDmerged_Total_aboveMedian.jpg"),
#       width = 18, height = 24, units = 'cm', dpi = 300)



ggsave(paste0(outpath,"2025_hapibd_matrixAdjustSize_length_sharing_IBDmerged_all.pdf"), width =24, height = 16, units = "cm",useDingbats=FALSE) # Figure 4A
ggsave(paste0(outpath,"2025_hapibd_matrixAdjustSize_length_sharing_IBDmerged_aboveMedian.pdf"), width =24, height = 16, units = "cm",useDingbats=FALSE) # Figure 4A




