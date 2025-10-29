
#### SCRIPT FOR GETTING MIXTURE MODEL and SOURCEFIND RESULTS
### Author: Victor Borda
### Created: April 2019

library(roxygen2)
library(xtable)
library(devtools)

sourcefind_path <-"/Users/vborda/Documents/Public_data/INS_data/Chromopainter_runs/contributions/sourcefind_20240904/"


### COMBINE MIXTURE MODEL RESULTS

pops <- regions<-c("Iquitos","Tumbes","Lambayeque","Trujillo","Ancash","PEL","Lima","Ica","Arequipa","Ayacucho","Moquegua","Tacna","Cusco","Puno")


##########################################
### AVERAGE AND COMBINE SOURCEFIND RESULTS
#sourcemat <- matrix(0,nc=length(pops),nr=length(mgpops))
#rownames(sourcemat) <- mgpops
#colnames(sourcemat) <- pops
mgpops <- pops[1:69] ###total number of population non aymaras, Shipibos, Quechua and two Chinese
#mgpops <- pops[1:72] ###total number of population to extract Con Quechuas and Aymaras
#mgpops <- pops[1:70] ###total number of population to extract
sourcemat <- matrix(0,nc=length(mgpops),nr=length(pops))
rownames(sourcemat) <- pops
colnames(sourcemat) <-c("Adygei", "Ashaninkas", "Awajun", "BEB", "Balochi", "BantuKenya", "Basque", "Bedouin",
                        "BergamoItalian", "Biaka", "Brahui", "Burusho", "CDX", "CEU", "CHB", "CHS", "Cambodian",
                        "Candoshi", "Chopccas", "Daur", "Druze", "ESN", "FIN", "French",
                        #"Candoshi", "Chopccas", "Daur", "Druze", "ESN", "FIN", "French",
                        "GBR", "GIH", "GWD", "Han", "Hazara", "Hezhen", "IBS", "ITU", "JPT", "KHV", "Kalash", "Karitiana",
                        "LWK", "Lahu", "Lamas", "MSL", "Makrani", "Matses", "Mbuti", "Miao", "Moche", "Mongolian", "Mozabite",
                        "Naxi", "Orcadian", "Oroqen", "PJL", "Palestinian", "Pathan", "Qeros", "Russian", "STU", "Sardinian",
                        "She", "Shimaa", "Sindhi", "Surui", "TSI", "Tallanes", "Tu", "Tujia", "Uros", "YRI", "Yakut",
                        "Yi")


## Results are average among rows for each column
for(pop in pops)
{
  #infile <- paste0(sourcefind_path,"/sourcefind_2023_non_QA_rephased_",pop,".txt")
  infile <- paste0(sourcefind_path,"/sourcefind_20240904_",pop,".txt")
  tablesource <- read.table(infile,header=T)
  sources_values<-tablesource[,3:69] ## Total number of columns with Quechuas and Aymaras
  #sources_values<-tablesource[,3:72] ## Total number of columns
  meansource<-colMeans(sources_values)
  sourcemat[pop,colnames(sources_values)] <- as.numeric(meansource)
  
}

#sourcemat <- read.table("/Users/vborda/Documents/Public_data/INS_data/Chromopainter_runs/contributions/Admixed_1.9M_sourcefind_2024.txt",header=T,row.names=1,as.is=T)


sourcemat[sourcemat<0.01] <- 0
sourcemat <- sourcemat/rowSums(sourcemat)

library(reshape2)
library(dplyr)

sourcemat.subset <- sourcemat[, colSums(sourcemat != 0) > 0]
sourcemat.subset <- tibble::rownames_to_column(as.data.frame(sourcemat.subset), "Pop")

sourcemat.subset$Pop = factor(sourcemat.subset$Pop, levels = c("Iquitos","Tumbes","Lambayeque","Trujillo","Ancash","PEL","Lima","Ica","Arequipa","Moquegua","Tacna","Ayacucho","Cusco","Puno"))


write.table(sourcemat.subset,paste0(sourcefind_path,"Admixed_1.9M_sourcefind_20250808.txt"), col.names=T,row.names=T)


mdat1 = melt(sourcemat.subset, id.vars="Pop", 
             variable.name="Ancestry", value.name="Proportion")

colors<-data.frame(

 ## non Quechuas and Aimaras
  #col1 = c("navy","dodgerblue3",
  #         "firebrick3","firebrick4","firebrick1",
  #         "yellow1","wheat1","yellow4",
  #         "palegreen", "green", "palegreen3","springgreen2","tan4","tan2","chocolate4", "lightskyblue3","powderblue" ),         
  
  #Ancestry =c("YRI","LWK",
  #            "IBS", "BergamoItalian","French",
  #            "CHB","Tu","Mongolian",
  #            "Awajun","Candoshi","Lamas", "Shipibo","Chopccas","Qeros", "Uros","Moche","Tallanes" )

  
 #with Quechuas and Aimaras
  # col1 = c("navy","dodgerblue",
  #           "firebrick1","firebrick3","red4","yellow1","yellow4","wheat1",
  #       "darkgreen","palegreen", "green", "palegreen3","springgreen2",
  #       "orange", "salmon2","chocolate4", "lightskyblue3","powderblue" ),         

   col1 = c("navy","dodgerblue",
             "darkorchid","mediumorchid","darkviolet","yellow1","yellow4","wheat1",
         "darkgreen","palegreen", "green", "palegreen3","springgreen2",
         "orange", "salmon2","chocolate4", "lightskyblue3","powderblue" ),
  
    Ancestry =c("YRI","LWK",
              "IBS","French","BergamoItalian", "CHB","Mongolian","Tu",
              "Ashaninkas","Awajun","Candoshi","Lamas","Matses",
              "Chopccas","Qeros", "Uros", "Moche","Tallanes" )

)

Ancestry =c("YRI","LWK",
            "IBS","French","BergamoItalian", "CHB","Mongolian","Tu",
            "Ashaninkas","Awajun","Candoshi","Lamas","Matses",
            "Chopccas","Qeros", "Uros", "Moche","Tallanes" )

#Ancestry =c("YRI","LWK",
#            "IBS", "BergamoItalian","French",
#            "CHB","Tu","Mongolian",
#            "Awajun","Candoshi","Lamas", "Shipibo","Chopccas","Qeros", "Uros","Moche","Tallanes" )

colors$Ancestry=factor(colors$Ancestry,levels = Ancestry)

mdat2<-merge(mdat1,colors, by.x = "Ancestry",by.y = "Ancestry")

mdat2$Ancestry=factor(mdat2$Ancestry,levels = Ancestry)

#mdat2<- mdat2 %>% arrange(Population) 

library(ggplot2)
library(forcats)
library(ggthemes)
library(patchwork)
outpath<-"/Users/vborda/Documents/Public_data/INS_data/Manuscript/R2_Figures"
  #jpeg(filename = paste0(outpath,"/Fig_2_C_20250808_SOURCEFIND.jpg"),width = 25, height = 15, units = "cm", pointsize =12,  res = 400)
  #jpeg(filename = paste0(sourcefind_path,"/SOURCEFIND_PGP_2023_ver_rephased_nonQA.jpg"),width = 25, height = 15, units = "cm", pointsize =12,  res = 400)
 sourcefind<-ggplot(mdat2, aes(Pop, Proportion, fill = Ancestry,order=Ancestry)) +
    geom_bar(stat="identity", position="fill", width=5, colour="grey25") +
    facet_grid(. ~ Pop, drop=TRUE, space="free", scales="free") +
    theme(panel.grid=element_blank()) +
    theme(panel.background=element_rect(fill=NA, colour=NA)) +
    theme(panel.margin.x=grid:::unit(0, "lines")) +
    theme(axis.title.x=element_blank()) +
    theme(axis.text.x=element_blank()) +
    theme(axis.ticks.x=element_blank()) +
    theme(strip.background=element_blank()) +
    theme(strip.text=element_text(size=12,angle = 90)) +
    theme(legend.key=element_rect(colour="grey25")) +
    scale_x_discrete(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0)) +
    scale_fill_manual(values=colors$col1) +
    
    scale_color_continuous(guide=guide_legend(direction = "right", title.position = "left",title.theme = element_text(angle = 90),
                                         label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))
    #guides(fill=guide_legend(override.aes = list(size = 1), direction = "horizontal", title.position = "left", title.theme = element_text(angle = 90),
    #                     label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
    
    #                     label.theme = element_text(angle = 90)))
     
  
  
  #dev.off()
  ggsave(paste0(outpath,"/Fig_1_C_20251027_SOURCEFIND.jpg"), sourcefind, bg="transparent",width = 25, height = 15, units = "cm")
  
  write.table(sourcemat,paste0(sourcefind_path,"Admixed_1.9M_sourcefind_20250808.txt"), col.names=T,row.names=T)


  
  
  
  
  
  

##########################################
### SAVING BOTH TABLES

#write.table(mixmat,paste0(results_path,"Admixed_1.9M_NNLS_May2019.txt"),col.names=T,row.names=T)



## PLOTTING BOTH RESULTS 

pcolshex <-popkey$Colour
poplabpos <- 0

pdf(paste0(results_path,"Admixed_1.9M_filtered_MM_SF.pdf"),height=5,width=7)
layout(matrix(c(1,2,3,4),byrow = TRUE, ncol = 4), #matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,6,7),3,5)
       widths=c(1,2,2,2))
topmar <- 4
n_pops <- 13
cores<-as.character(pcolshex)
bp0 <- barplot(t(as.matrix(x1)),
              yaxt="n",xlab="",beside=F,
              main="",cex.main=0.75,border=NA,horiz=T,
              col=cores,xaxt="n",add=F,plot=F) ## MODIFICADO col=srccols
 #par(mar=c(6,0,topmar,0.5)) AQUI SE MODIFICA EL GROSOR DE LOS PLOTS
par(mar=c(1.5,0,topmar+0.5,0)) 
plot(0,0,xlim=c(0,1), ylim=c(bp0[length(bp0)],bp0[1]),type="n",axes=F,xlab="",ylab="")

#text(y=-0.9,x=-0.015,labels="Sources involved in admixture",xpd=T,srt=0,adj=0,cex=1)  #
#text(y=-1,x=0.2,labels="MIXTURE\n  MODEL",xpd=T,srt=0,adj=0)

text(1,0.6,labels="Tumbes", col="cornflowerblue",las=2,tck=0,cex=1.1,pos=2)
text(1,1.8,labels="Lambayeque", col="cornflowerblue",las=2,tck=0,cex=1.1,pos=2)
text(1,3.1,labels="Ancash", col="darksalmon",las=2,tck=0,cex=1.1,pos=2)
text(1,4.3,labels="Trujillo", col="cornflowerblue",las=2,tck=0,cex=1.1,pos=2)
text(1,5.5,labels="Lima", col="cornflowerblue",las=2,tck=0,cex=1.1,pos=2)
text(1,6.7,labels="Afrodescendants", col="cornflowerblue",las=2,tck=0,cex=1,pos=2)
text(1,7.9,labels="Moquegua", col="cornflowerblue",las=2,tck=0,cex=1.1,pos=2)
text(1,9.1,labels="Tacna", col="cornflowerblue",las=2,tck=0,cex=1.1,pos=2)
text(1,10.2,labels="Arequipa", col="darksalmon",las=2,tck=0,cex=1.1,pos=2)
text(1,11.5,labels="Ayacucho", col="darksalmon",las=2,tck=0,cex=1.1,pos=2)
text(1,12.7,labels="Cusco", col="darksalmon",las=2,tck=0,cex=1.1,pos=2)
text(1,14,labels="Puno", col="darksalmon",las=2,tck=0,cex=1.1,pos=2)
text(1,15.2,labels="Iquitos", col="olivedrab4",las=2,tck=0,cex=1.1,pos=2)

x1 <- mixmat[pops[1:13],pops]
bp <- barplot(t(as.matrix(x1)),
              yaxt="n",xlab="",beside=F,
              main="",cex.main=0.75,border=NA,horiz=T,
              col=cores,xaxt="n",add=F,plot=F) ## MODIFICADO col=srccols
par(mar=c(1.5,0.5,topmar+0.5,0)) #par(mar=c(6,0,topmar,0.5)) AQUI SE MODIFICA EL GROSOR DE LOS PLOTS
plot(0,0,xlim=c(0,1), ylim=c(bp[length(bp)],bp[1]),type="n",axes=F,xlab="",ylab="")

bp <- barplot(t(as.matrix(x1)),width = 1,
              yaxt="n",xlab="",beside=F,
              main="",cex.main=0.75,border=NA,horiz=T,
              col=cores,xaxt="n",add=T,plot=T)

#text(y=-0.9,x=-0.015,labels="Sources involved in admixture",xpd=T,srt=0,adj=0,cex=1)  #
text(y=-1,x=0.25,labels="MIXTURE MODEL",xpd=T,srt=0,adj=0,cex = 1.2)




x2 <- sourcemat[pops[1:13],pops]
bp2 <- barplot(t(as.matrix(x2)),
              yaxt="n",xlab="",beside=F,
              main="",cex.main=0.75,border=NA,horiz=T,
              col=cores,xaxt="n",add=F,plot=F) ## MODIFICADO col=srccols
par(mar=c(1.5,0.5,topmar+0.5,0)) #par(mar=c(6,0,topmar,0.5)) AQUI SE MODIFICA EL GROSOR DE LOS PLOTS
plot(0,0,xlim=c(0,1), ylim=c(bp2[length(bp2)],bp[1]),type="n",axes=F,xlab="",ylab="")

bp2 <- barplot(t(as.matrix(x2)),width = 1,
              yaxt="n",xlab="",beside=F,
              main="",cex.main=0.75,border=NA,horiz=T,
              col=cores,xaxt="n",add=T,plot=T)
text(y=-1,x=0.25,labels="SOURCEFIND",xpd=T,srt=0,adj=0,cex=1.2)
## LEGEND
par(mar=c(1,0,4,0))
plot(0,0,axes=F,xlab="",ylab="",type="n")
legend_text <- c("Tallanes","Moche","Jacarus","Quechuas","Chopccas","Qeros","Aimaras","Uros","Ashaninkas","Matsiguenkas","Matses","Nahua","Shipibo","Candoshi","Awajun","Lamas","Chachapoyas","South Europe","North Europe","West Africa","West Central Africa","East Africa","East Asian")
l <- legend("top",legend=legend_text, pch=c(rep(15,10)),
            col=c("cornflowerblue","cadetblue","peachpuff4","tan","tan1","tan2","tan3","tan4","darkolivegreen","darkolivegreen1","darkolivegreen2","darkolivegreen3","chartreuse","chartreuse3","darkolivegreen4","greenyellow","slategray2",
"firebrick1","firebrick3","blue","blue2","blue4","gold"),bty="n",
            ncol=1,xpd=T,pt.cex=1,x.intersp=1,y.intersp=1,
            pt.lwd=2,cex=1, title="ANCESTRY REGION")


###################################################################

dev.off()





