### plotting HApNE

path<-"/Users/vborda/Documents/Public_data/INS_data/IBD_Analyses/HapNe/"
regions<-c("Iquitos","Tumbes","Lambayeque","Trujillo","Huaraz","PEL","Lima","El Carmen","Arequipa","Ayacucho","Moquegua","Tacna","Cusco","Juliaca")

library(stringr)
library(dplyr)
library(gapminder)
library(gghalves)
library(ggdist)
library(ggpubr)
library(ggforce)
library(reshape2)
library(ggplot2)
library(ggbreak)
library(patchwork)

for (Pop in c(1:length(regions))){
  Ne <- read.table(file=paste0(path,"202408_Results/",regions[Pop],"/hapne.csv"),header=TRUE,sep=",")
  #toplot<-Ne[6:51,1:4]
  
  p<-ggplot(Ne, aes(x=TIME, y=Q0.5)) +
    geom_ribbon(aes(ymin = Q0.025, ymax = Q0.975), alpha=0.2, fill="navy") +
    geom_ribbon(aes(ymin = Q0.25, ymax = Q0.75), alpha=0.2,fill="dodgerblue") +
    geom_line(col="firebrick") +
    ggtitle(regions[Pop]) +
    scale_y_log10(limits = c(9e2,6e7)) +
     xlab("Generations ago") + ylab("Effective Population Size")  +
    theme(plot.title = element_text(hjust = 0,size = 8, face = "bold"), axis.text = element_text(size = 6),
          axis.title=element_text(size=8),plot.margin = unit(c(0.7,0.7,0.7,0.7), 'lines')) #+
    #scale_y_continuous(trans='log10') 
    #coord_cartesian(ylim = c(0, 500000)) 
  #+ scale_y_break(c(120000, 1500000)) +
    #scale_y_continuous(breaks=c(0,50000, 100000,150000,1500000,1600000),limits =c(0,1600000) )
  
  nam <- paste("Ne_", regions[Pop], sep = "")
  assign(nam, p)
  
}

library(gridExtra)

figure <- ggarrange(Ne_Tumbes,Ne_Lambayeque,Ne_Iquitos,Ne_Trujillo,Ne_Ancash,Ne_Lima,
                       ncol = 2, nrow = 3)
jpeg(filename = paste0(path,"Test_Ne.jpg"),width = 18, height = 18, units = "cm", pointsize =12,  res = 300) 
figure
dev.off()


figure2 <- ggarrange(Ne_Ica,Ne_Arequipa,Ne_Ayacucho,Ne_Cusco,Ne_Moquegua,Ne_Tacna,Ne_Puno,
                    ncol = 3, nrow = 3)
jpeg(filename = paste0(path,"Test_Ne2.jpg"),width = 18, height = 18, units = "cm", pointsize =12,  res = 300) 
figure2
dev.off()


full<-ggarrange(Ne_Iquitos,Ne_Tumbes,Ne_Lambayeque,Ne_Trujillo,Ne_Huaraz,Ne_PEL,Ne_Lima,
                `Ne_El Carmen`,Ne_Ayacucho,Ne_Arequipa,Ne_Moquegua,Ne_Tacna,Ne_Cusco,Ne_Juliaca,
                ncol = 4, nrow = 4)
outpath<-"/Users/vborda/Documents/Public_data/INS_data/Manuscript/R1_Figures/"
jpeg(filename = paste0(outpath,"Fig_S9_202508_HapNe_LOG.jpg"),width = 21, height = 18, units = "cm", pointsize =12,  res = 300) 
full
dev.off()


### no LOG


for (Pop in c(1:length(regions))){
  #Pop<-3
  Ne <- read.table(file=paste0(path,"202408_Results/",regions[Pop],"/hapne.csv"),header=TRUE,sep=",")
  #toplot<-Ne[6:51,1:4]
  
  p<-ggplot(Ne, aes(x=TIME, y=Q0.5)) +
    geom_ribbon(aes(ymin = Q0.025, ymax = Q0.975), alpha=0.2, fill="navy") +
    geom_ribbon(aes(ymin = Q0.25, ymax = Q0.75), alpha=0.2,fill="dodgerblue") +
    geom_line(col="firebrick") +
    ggtitle(regions[Pop]) +
    scale_y_break(c(120000, 1500000)) +
    #scale_y_log10(limits = c(9e2,3e6)) +
    xlab("Generations ago") + ylab("Effective Population Size")  +
    theme(plot.title = element_text(hjust = 0,size = 8, face = "bold"), axis.text = element_text(size = 6),
          axis.title=element_text(size=8),plot.margin = unit(c(0.7,0.7,0.7,0.7), 'lines')) #+
  #scale_y_continuous(trans='log10') 
  #coord_cartesian(ylim = c(0, 500000)) 
  #+ scale_y_break(c(120000, 1500000)) +
  #scale_y_continuous(breaks=c(0,50000, 100000,150000,1500000,1600000),limits =c(0,1600000) )
  
  nam <- paste("Ne_", regions[Pop], sep = "")
  assign(nam, p)
  
}


full<-ggarrange(Ne_Iquitos,Ne_Tumbes,Ne_Lambayeque,Ne_Trujillo,Ne_Huaraz,Ne_PEL,Ne_Lima,
                `Ne_El Carmen`,Ne_Ayacucho,Ne_Arequipa,Ne_Moquegua,Ne_Tacna,Ne_Cusco,Ne_Juliaca,
                ncol = 4, nrow = 4)

#jpeg(filename = paste0(path,"2025January_HapNe_NoLOG.jpg"),width = 21, height = 18, units = "cm", pointsize =12,  res = 300) 
#full
#dev.off()

outpath<-"/Users/vborda/Documents/Public_data/INS_data/Manuscript/R1_Figures/"
jpeg(filename = paste0(outpath,"Fig_S10_202508_HapNe_NoLOG.jpg"),width = 21, height = 18, units = "cm", pointsize =12,  res = 300) 
full
dev.off()



##########################
#### unique plot
##########################

path<-"/Users/vborda/Documents/Public_data/INS_data/IBD_Analyses/HapNe/"
#regions<-c("Tumbes","Lambayeque","Iquitos","Puno","Lima","Chincha","Trujillo","Huaraz","Arequipa","Huamanga","Mariscal Nieto","Cusco","Tacna")

multiple.Ne <- matrix(0,nc=7,nr=0)
colnames(multiple.Ne)<-c("TIME","Q0.025","Q0.25","Q0.5","Q0.75","Q0.975","Population")

for (Pop in c(1:length(regions))){
  Ne <- read.table(file=paste0(path,"202408_Results/",regions[Pop],"/hapne.csv"),header=TRUE,sep=",")
  
  Ne$Population<-regions[Pop]
  
  multiple.Ne<-rbind(multiple.Ne,Ne)
  
}

multiple.Ne$Population <- factor(multiple.Ne$Population,
                                 levels = c("Iquitos","Tumbes","Lambayeque","Trujillo","Huaraz","PEL","Lima","El Carmen",
                                            "Arequipa","Ayacucho","Moquegua","Tacna","Cusco","Juliaca"))

multiple.Ne$Linetype <- ifelse(as.numeric(as.factor(multiple.Ne$Population)) %% 2 == 0, "solid", "dashed")


p<-ggplot(multiple.Ne, aes(TIME, Q0.5, group = Population)) + 
  geom_ribbon(aes(ymin = Q0.025, ymax = Q0.975, fill = Population), alpha = 0.3) + 
  geom_line(aes(color = Population, linetype = Population), alpha = 0.6) + 
  scale_y_log10(limits = c(7e2, 2e8)) + 
  xlim(0, 40) + 
  xlab("Generations ago") + 
  ylab("Log Effective Population Size") + 
  scale_color_manual(values = scales::hue_pal()(length(unique(multiple.Ne$Population)))) + 
  scale_fill_manual(values = scales::hue_pal()(length(unique(multiple.Ne$Population)))) + 
  scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), length.out = length(unique(multiple.Ne$Population)))) + 
  guides(
    color = guide_legend(title = "Population", override.aes = list(size = 2)),  # Increase line thickness in legend
    linetype = guide_legend(title = "Population", override.aes = list(size = 2)) # Increase line thickness in legend
  ) + 
  theme(
    legend.title = element_text(size = 8), # Adjust legend title size
    legend.text = element_text(size = 7),  # Adjust legend item text size
    legend.key.size = unit(0.5, "cm"),     # Adjust size of legend keys
    legend.spacing.y = unit(0.1, "cm"),    # Reduce vertical spacing between legend items
    plot.margin = unit(c(1,3,1,3), 'lines') #bottom, left, top, right
  )

outpath
jpeg(filename = paste0(outpath,"2025January_HapNe_allpops.jpg"),width = 20, height = 12, units = "cm", pointsize =12,  res = 300) 
p
dev.off()

load(file="/Users/vborda/Documents/Public_data/INS_data/Manuscript/R2_Figures/SubFigures/Figure_3_above_ggMedian_8G.RData")

above_row_all <- ggarrange(gg_older, gg_8G, ncol=2, labels = c("A","B"), vjust = 2.2,
                        font.label=list(color="black",size=11))


above_row <- ggarrange(gg_median_older, gg_median_8G, ncol=2,nrow=1, labels = c("A","B"), vjust = 2.2,
                       font.label=list(color="black",size=11))


full<-ggarrange(above_row,p,
                labels = c("", "C"),heights=c(0.8,1),
                ncol = 1, nrow = 2, font.label=list(size=11))

output_path<-"/Users/vborda/Documents/Public_data/INS_data/Manuscript/R1_Figures/"
jpeg(filename = paste0(output_path,"Figure_3_202508_IBDSharing_HapNe.jpg"),width = 20, height = 18, units = "cm", pointsize =12,  res = 600) 
full
dev.off()

output_path<-"/Users/vborda/Documents/Public_data/INS_data/Manuscript/R1_Figures/"
jpeg(filename = paste0(output_path,"Figure_3_EXTENDED_202508_IBDSharing_NewSup.jpg"),width = 22, height = 8, units = "cm", pointsize =12,  res = 300) 
above_row_all
dev.off()


path<-"/Users/vborda/Documents/Public_data/INS_data/Manuscript/R2_Figures/"
tiff(filename = paste0(path, "Figure_3_202510_IBDSharing_HapNe.tiff"),
     width = 20, height = 18, units = "cm",
     res = 300, pointsize = 12)

full

dev.off()
