
# This pipeline runs R's MDS on the masked data.
# Make sure the masked (missing) data are coded as NA's.
# We have run this approach with over 10,000 individuals and 10,000 SNPs: to do so took around 1 day of computing time.
# Computing with fewer individuals is faster. 
path<-"/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/"
for(ancestry in 2){
  gts=read.table(paste0(path,"allchr_thinned_anc",ancestry,"_apropgt0.5_singleref.bgl"),header=T)
  #rownames(gts)=gts[,2]
  rownames(gts)=paste0("rs",seq(1,nrow(gts)))
  gts=gts[,-c(1,2)]
  gts=t(gts)
  # Remove low frequency variants.
  maf = apply(gts,2,function(z){
    vars=unique(z)
    vars=vars[!is.na(vars)]
    if(length(vars)!=2)return(0)
    n1 = sum(z==vars[1],na.rm=T)
    n2 = sum(z==vars[2],na.rm=T)
    return(min(n1,n2)/(n1+n2))
  })
  gts=gts[,maf>0.1]
  finaldata=apply(gts,2,function(zz)as.numeric(as.factor(zz)))
  rownames(finaldata)=rownames(gts)
  ## At this point you may want to remove some of the individuals before proceeding to reduce computing time and to focus the PCs on the individuals of interest.
  finaldata=data.frame(finaldata)
  # Calculate distance matrix.
  mydist=dist(finaldata)
  # Fix NA's in distance matrix that are due to some haplotypes not having overlapping positions without masking.
  mydist2=mydist
  mydist2[is.na(mydist2)]=mean(mydist,na.rm=T)
  # Calculate 4 PCs (other numbers of PCs are possible).
  mymds=cmdscale(mydist2,10)
  # Output the results to a text file.
  write.table(mymds,paste("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/Ancestry",ancestry,"_specific_pcs",sep=""), quote=F,col.names=F)
}

## The PCs written out to the file can be used to create plots with R or with another package.

#pca<-read.table("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/Ancestry2_specific_pcs", header=F)
pca<-read.table("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/Ancestry4_specific_pcs", header=F) ## Natives

colnames(pca)<-c("ID","PC1","PC2","PC3","PC4")
path<-"/Users/vborda/Documents/Public_data/INS_data/data+RefPops/"
demo<- read.table(paste0(path,"/FULL_1KGP_HGDP_PGP_IDs.txt"), head=FALSE, fill = TRUE,comment.char = "")
key<- read.table("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/pop_colors_EUR.txt", head=TRUE, fill = TRUE,comment.char = "")
##colnames(demo)<-c("ID","POP")
colnames(demo)<-c("ID","POP","tags","labels")
pca$ID2 <- sub('.1$','',pca$ID)
toplot<-merge(pca,demo, by.x = "ID2", by.y = "ID")
toplot2<-merge(toplot,key, by.x = "POP", by.y = "POP")

library(ggplot2)
library(gapminder)
library(gghalves)
library(ggdist)
library(ggpubr)
library(ggforce)
library(ggmagnify)
library(ggfx)


regions<-data.frame(
  Deparment=c("Tumbes","Lambayeque","Iquitos","Trujillo","Ancash","Lima","Ica","Arequipa","Ayacucho","Moquegua","Tacna","Cusco","Puno","PEL"),
  Province=c("Tumbes","Lambayeque","Iquitos","Trujillo","Huaraz","Lima","Chincha","Arequipa","Huamanga","Mariscal Nieto","Tacna",
             "Cusco / Paucartambo","Puno/San Roman","PEL"),
  tag=c(1:14),
  symbol=c(7,21,8,22,9,23,10,24,12,25,13,4,14,15),
  color=c("cornflowerblue","darkred","darkolivegreen1","darkorange","cyan","purple","gold1",
          "cornflowerblue","darkred","darkolivegreen1","darkorange","cyan","purple","tan")
  )

REF<-data.frame(
  Province=c("Tallanes","Moche","Chachapoyas","Awajun","Candoshi","Matses","Lamas","Ashaninkas","Matsiguenkas","Shimaa","Nahua","Colombian",
             "Surui","Karitiana","Shipibo",
             "Aimaras","Chopccas","Jacarus","Qeros","Quechuas","Uros",
             "Pima","Maya","MXL","CLM"),
  Area=c("Coast","Coast","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian",
         "Amazonian","Amazonian","Amazonian",
         "Andean","Andean","Andean","Andean","Andean","Andean",
         "MesoAmerican","MesoAmerican","Admixed","Admixed"),
  symbol=c(15,15,rep(16, times = 13),rep(17, times = 6),18,18,6,6)
)
REF$color<-"grey"

#remove<-c("Tumbes","Lambayeque","Iquitos","Puno","Lima","Ica","Trujillo","Ancash","Arequipa","Ayacucho","Moquegua","Cusco","Tacna",
#           "CLM","PEL","PUR","Shipibo","MXL")

references<-toplot2[!(toplot2$POP %in% regions$Deparment),]
references$alpha<-0.2
#references$symbol<-19

targets<-toplot2[(toplot2$POP %in% regions$Deparment),]
targets$alpha<-1
targets<-merge(targets[,c(1:9)],regions,by.x="POP",by.y="Deparment")
references<-merge(references[,c(1:9)],REF,by.x="POP",by.y="Province")

library(tidyverse)
library(palmerpenguins)



arrows

p1<-ggplot() +              
  geom_point(data=references,aes(PC1,-PC2),shape=references$symbol,alpha=0.5,color=references$color,size=3) +
  geom_point(data=targets,aes(PC1,-PC2),shape=targets$symbol,alpha=1,color=targets$color,size=2)+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  ) 

 theme(
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.background = element_rect(fill='transparent'),
         legend.box.background = element_rect(fill='transparent')
       )


arrows <- 
  tibble( #PEL
    x1 = c(-18,30, 26,-7.5,-11,-15,-17,-16,3,-1,-4,0.8,-8,-10,-20,35,35),
    x2 = c(-20,21, 4,0,-4,-9.5,-9.5,-12,-10,-5,-12,-14,-12.5,-16,-16,40,28),  # HEAD
    y1 = c(-1,8.5, 19.5,24,20,16,10,5.7,-2,-6.1,-10.8,-16,-20,-25,-29,-29,-29),
    y2 = c(0,4, 15,17.5,10,6,3.5,3,-5.5,-6.5,-9.5,-14,-16,-17.5,-20.5,-16,-2)
  )

p2<-p1 +  
  annotate("text", x = 30, y = 10, label = "Iquitos",size=5) +
  annotate("text", x = 27, y = 21, label = "Tumbes",size=5)+
  annotate("text", x = -15, y = 25, label = "Lambayeque",size=5)+
  annotate("text", x = -15, y = 21, label = "Trujillo",size=5)+
  annotate("text", x = -15, y = 17, label = "Huaraz",size=5)+
  annotate("text", x = -19, y = 11.5, label = "Chincha",size=5)+
  annotate("text", x = -20, y = 6, label = "Lima",size=5)+
  annotate("text", x = 5, y = -1, label = "Huamanga",size=5)+
  annotate("text", x = 6, y = -6, label = "Arequipa",size=5)+
  annotate("text", x = 8, y = -11, label = "Mariscal Nieto",size=5)+
  annotate("text", x = 6, y = -16, label = "Tacna",size=5)+
  annotate("text", x = 5, y = -21, label = "Puno / San Roman",size=5)+
  annotate("text", x = 3, y = -26, label = "Cusco / Paucartambo",size=5)+
  annotate("text", x = 35, y = -30.5, label = "AMAZONIAN", colour = "darkgreen", fontface =2,size=5)+
  annotate("text", x = -17, y = -30.5, label = "ANDEAN", colour = "tan4", fontface =2,size=5)+
  annotate("text", x = 2, y = 26, label = "COAST", colour = "dodgerblue", fontface =2,size=5)+
  geom_curve(
    data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.08, "inch"),type="closed"), size = 0.5,
    color = "gray20", curvature = 0)


p2<-p1+ geom_curve(
  data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2),
  arrow = arrow(length = unit(0.08, "inch"),type="closed"), size = 0.5,
  color = "gray20", curvature = 0) + theme(panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent'))
    


jpeg("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/20241031_NAT_PCA_nolabes.jpg",width = 15, height = 15, units = "cm", pointsize =12,  res = 600) 
p2
dev.off()

##### EDITING PLOT 20250808
# Separate PEL from other targets
pel_points <- targets %>% filter(POP == "PEL")
other_targets <- targets %>% filter(POP != "PEL")

# Plot order: PEL first (behind), then others, then references
p1 <- ggplot() +              
  geom_point(data = references,
             aes(PC1, -PC2),
             shape = references$symbol,
             alpha = 0.5,
             color = references$color,
             size = 3) +
  
  geom_point(data = pel_points,  # PEL in background
             aes(PC1, -PC2),
             shape = pel_points$symbol,
             alpha = 1,
             color = pel_points$color,
             size = 2) +
  
  geom_point(data = other_targets,  # others on top
             aes(PC1, -PC2),
             shape = other_targets$symbol,
             alpha = 1,
             color = other_targets$color,
             size = 2) +
  
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )

arrows <- 
  tibble( #PEL,PEL
    x1 = c(-20,-20,30, 26,-7.5,-11,-15,-17,-16,3,-1,-4,0.8,-8,-10,-20,35,35),
    x2 = c(-14.5,-18,21, 4,0,-4,-9.5,-9.5,-12,-10,-5,-12,-14,-12.5,-16,-16,40,28),  # HEAD
    y1 = c(3,0,8.5, 19.5,24,20,16,10,5.7,-2,-6.1,-10.8,-16,-20,-25,-29,-29,-29),
    y2 = c(4.5,-1,4, 15,17.5,10,6,3.5,3,-5.5,-6.5,-9.5,-14,-16,-17.5,-20.5,-16,-2)
  )

p2<-p1+ geom_curve(
  data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2),
  arrow = arrow(length = unit(0.08, "inch"),type="closed"), size = 0.5,
  color = "gray20", curvature = 0) + theme(panel.background = element_rect(fill='transparent'),
                                           plot.background = element_rect(fill='transparent', color=NA),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           legend.background = element_rect(fill='transparent'),
                                           legend.box.background = element_rect(fill='transparent'))

ggsave("/Users/vborda/Documents/Public_data/INS_data/Manuscript/R1_Figures/Fig_2_B_20250808_asPCA.png", p2, bg="transparent",width = 15, height = 15, units = "cm")
#ggsave("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/20241031_NAT_PCA_nolabes.png", p2, bg="transparent",width = 15, height = 15, units = "cm")















#targets$symbol<-15
#targets$color<-"grey"

#targets2<-merge(targets,regions,by.x="POP",by.y="Deparment")
toplot3<-rbind(references,targets)

#targets3<-merge(toplot2,regions,by.x="POP",by.y="Deparment",all.x=TRUE)

targets3<-merge(toplot3,regions,by.x="POP",by.y="Deparment",all.x=TRUE)



#tiff(paste0(path,"/Local_ancestry/aspca/PCA_AFR.tiff"), width = 15, height = 15, units = 'in', res = 300)
#ggplot(toplot, aes(x=PC1, y=PC2,color=POP))  + geom_text(size=3,label=toplot$POP) + theme(legend.position = "none") + #+ geom_point(shape=2 )
#  geom_magnify(from = c(-10, 12, -12, 10), 
#               to = c(30, 120, 10, 90), 
#               shadow = TRUE)

#dev.off()

### EUR
references$POP <- factor(references$POP, levels=c('FIN', 'Russian', 'GBR', 'Orcadian','CEU','IBS','French','Basque',
                                                   'TSI','Tuscan','BergamoItalian','Sardinian','Adygei','MXL', 'PUR', 
                                                   'CLM','PEL','Shipibo'))
scalesymbols <- c(24, 17, 22, 22,22,19,25,13,
                  8,19,19,23,19,18,19,
                  19,19,19,
                  32,32,32,32,32,32,32,32)
scalecolors <- c('firebrick3', 'orange2', 'firebrick3', 'firebrick3','tomato1','magenta1','firebrick3','firebrick4',
                                                      'firebrick3','firebrick3','firebrick3','purple','purple','wheat', 'cyan', 
                                                      'grey61','grey61','green2',"black","black","black","black",
                                                      "black","black","black","black")

targets$POP <- as.factor(targets$POP)

colors=rbind(unique(references[,c("POP","color","symbol")]),unique(targets[,c("POP","color","symbol")]))

tiff("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/aspca_EUR_ver3.tiff", width = 20, height = 10, units = 'cm', res = 300)
ggplot(aes(PC1, PC2,color=as.factor(POP),shape=as.factor(POP)),data=references) +              
  geom_point(size=1.5,alpha=0.6,stroke=0.7) +
  geom_text(aes(PC1, PC2),data = targets,size=1.8,label=targets$POP) +
  scale_shape_manual(name="Populations",values = scalesymbols) +
  scale_color_manual(name="Populations",values = scalecolors) +
    facet_zoom(xlim=c(-15,0),ylim = c(-12,10),zoom.size = 1.2)+ theme_light() +
  theme(legend.text = element_text(size=5),
      legend.title=element_text(size=5),
      legend.position = "right",legend.key.size = unit(0.3, 'cm'))
dev.off()

tiff("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/aspca_EUR_ver3_pc3pc4.tiff", width = 20, height = 10, units = 'cm', res = 300)
ggplot(aes(PC3, PC4,color=as.factor(POP),shape=as.factor(POP)),data=references) +              
  geom_point(size=1.5,alpha=0.6,stroke=0.7) +
  geom_text(aes(PC1, PC2),data = targets,size=1.8,label=targets$POP) +
  scale_shape_manual(name="Populations",values = scalesymbols) +
  scale_color_manual(name="Populations",values = scalecolors) +
  facet_zoom(xlim=c(-15,0),ylim = c(-12,10),zoom.size = 1.2)+ theme_light() +
  theme(legend.text = element_text(size=5),
        legend.title=element_text(size=5),
        legend.position = "right",legend.key.size = unit(0.3, 'cm'))
dev.off()







### Natives
#p1<-ggplot(references, aes(PC1, -PC2)) +         
p1<-ggplot(targets3, aes(PC1, -PC2)) +         
  geom_point(color=targets3$color,shape=targets3$symbol,alpha=targets3$alpha,size=4) +
  #geom_point(aes(fill=color),shape=21,alpha=0.5,size=4) +
  geom_text(aes(label = tag), vjust = 0.5, hjust = 0.5, size=4) +
  coord_cartesian(xlim = c(-150, 60),ylim = c(-70, 25)) + 
  geom_magnify(from = c(-20, 10, -18, 20), 
               to = c(-150, -25, -70, 25), 
               shadow = TRUE,expand = FALSE) 
  #coord_cartesian(xlim = c(-40, 60),ylim = c(-70, 25)) + 
  #geom_magnify(from = c(-20, 0, -18, 10), 
  #             to = c(-40, 44, -70, -21), 
  #             shadow = TRUE,expand = FALSE) 

  #coord_cartesian(xlim = c(-55, 60),ylim = c(-180, 30)) + 
  #geom_magnify(from = c(-20, 25, -18, 20), 
  #             to = c(-55, 48, -180, -28), 
  #             shadow = TRUE,expand = FALSE) 





ggplot(references, aes(PC3, PC4)) +              
  geom_point(shape=references$symbol,alpha=0.3,color=references$color) +
  #geom_text(data = targets,size=2,label=targets$POP) +
  geom_text(data = targets,size=1.5,label=targets$tags)

#p1<-ggplot(references, aes(PC1, -PC2,shape=references$symbol,fill=references$color)) +              




ggplot() + 
  geom_point(data=df1, aes(x=day, y=sales), color='steelblue') + 
  geom_point(data=df2, aes(x=day, y=sales), color='coral2')


p2<-ggplot(references, aes(PC1, -PC2)) +              
  geom_point(shape=references$symbol,alpha=references$alpha,color=references$color) +
  geom_text(data = targets,size=0.7,label=targets$POP)+
  theme(
    panel.background = element_rect(fill='grey'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

p3<-ggplot(references, aes(PC1, -PC2)) +              
  geom_point(shape=references$symbol,alpha=references$alpha,color=references$color) +
  #geom_text(data = targets,size=1,label=targets$POP)+
  theme(
    panel.background = element_rect(fill='grey'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )



#tiff("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/202408_aspca_NAT_ver2.tiff", width = 15, height = 30, units = 'cm', res = 300)
tiff("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/202408_aspca_NAT_ver3.tiff", width = 15, height = 15, units = 'cm', res = 300)

  p1 #+ annotation_custom(ggplotGrob(p2), xmin = 32, xmax = 62, 
    #                     ymin = 0, ymax = 25)
dev.off()


tiff("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/aspca_NAT_ver4.tiff", width = 22, height = 18, units = 'cm', res = 300)

p1 + annotation_custom(ggplotGrob(p3), xmin = 32, xmax = 65, 
                       ymin = -1, ymax = 27)
dev.off()


tiff("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/aspca_NAT_ver3_small.tiff", width = 22, height = 18, units = 'cm', res = 150)

p1 + annotation_custom(ggplotGrob(p2), xmin = 35, xmax = 60, 
                       ymin = 0, ymax = 25)
dev.off()

### EAS
tiff("/Users/vborda/Documents/Public_data/INS_data/data+RefPops/Local_ancestry/aspca/aspca_EAS_ver2.tiff", width = 22, height = 22, units = 'cm', res = 150)

ggplot(references, aes(-PC2, PC1)) +              
  geom_point(shape=references$symbol,alpha=0.3,color=references$color) +
  geom_text(data = targets,size=2,label=targets$POP) 
dev.off()


library(remotes)
#install_github("r-spatial/sf")

library(sp)
library(maps)
library(maptools)
library(ggplot2)
library(ggthemes)

library(ggimage)
library(rnaturalearth)
library(countrycode)
library(ggthemes)
library(reshape2)

map <- ne_countries(scale = "medium", returnclass = "sf")
EAS<-map[map$subregion=="Eastern Asia",]
Cambodia<-map[map$admin=="Cambodia",]
Vietnam<-map[map$admin=="Vietnam",]


library(rvest)


coord <- read_html("https://developers.google.com/public-data/docs/canonical/countries_csv")

coord_tables <- coord %>% html_table(header = TRUE, fill = TRUE)

coord <- coord_tables[[1]]

