
library("GENESIS")
library(gdsfmt)
library(SNPRelate)
library(SeqVarTools)
library(GWASTools)
source("/Users/vborda/Documents/Public_data/INS_data/Phenotype_Associations/topmed.R")


## Running PCA

name.file<-"20231226_PGP_4Continental_groups_MAF_LD0.4"
path.samples<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/Inputs//"

########################

bed.fn <- paste0(path.samples,name.file,".bed")
bim.fn <- paste0(path.samples,name.file,".bim")
fam.fn <- paste0(path.samples,name.file,".fam")

snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, paste0(path.samples,"2024_GenomeWide_INS_LDGH_1KGP_HGDP.gds"))
snpgdsSummary(paste0(path.samples,"2024_GenomeWide_INS_LDGH_1KGP_HGDP.gds"))


### This data includes the whole high coverage 1KGP (Unrelated + Trios)

(gdsfile <- snpgdsOpen(paste0(path.samples,"2024_GenomeWide_INS_LDGH_1KGP_HGDP.gds")))
#snpgdsClose(gdsfile)

### LD Prunning

snpset <- snpgdsLDpruning(gdsfile, method="corr", slide.max.bp=10e6, 
                          ld.threshold=sqrt(0.4), verbose=FALSE)
pruned <- unlist(snpset, use.names=FALSE)
length(pruned)

### Pairwise Measures of Ancestry Divergence

ibd.robust  <- snpgdsIBDKING(gdsfile, num.thread=8)

KINGmat<-kingToMatrix(ibd.robust)

mypcair <- pcair(gdsfile, kinobj = KINGmat, divobj = KINGmat, kin.thresh=2^(-9/2), div.thresh=-2^(-9/2),
                 snp.include = pruned, maf=0.05)

#Identifying relatives for each sample using kinship threshold 0.0441941738241592
#Identifying pairs of divergent samples using divergence threshold -0.0441941738241592

#mypcair <- pcair(gdsfile, kinobj = KINGmat, divobj = KINGmat,
#                 snp.include = pruned)

save(mypcair,file=paste0(path.samples,"20240223_INS_LDGH_1KGP_HGDP_subset_PCA.RData"))

path.samples<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/Inputs//"
load(file=paste0(path.samples,"20240223_INS_LDGH_1KGP_HGDP_subset_PCA.RData"))


pca<-data.frame(ID = mypcair$sample.id,PC1=mypcair$vectors[,1],PC2=mypcair$vectors[,2],PC3=mypcair$vectors[,3],PC4=mypcair$vectors[,4],PC5=mypcair$vectors[,5],
                PC6=mypcair$vectors[,6],PC7=mypcair$vectors[,7],PC8=mypcair$vectors[,8],PC9=mypcair$vectors[,9],PC10=mypcair$vectors[,10],
                PC11=mypcair$vectors[,11],PC12=mypcair$vectors[,12],PC13=mypcair$vectors[,13],PC14=mypcair$vectors[,14],PC15=mypcair$vectors[,15],
                PC16=mypcair$vectors[,16],PC17=mypcair$vectors[,17],PC18=mypcair$vectors[,18],PC19=mypcair$vectors[,19],PC20=mypcair$vectors[,20])


#write.table(pca,"/Users/vborda/Documents/Public_data/INS_data/PCA_UMAP/20240219_PCA_INS_LDGH_1KGP_HGDP.txt", sep='\t', col.names = TRUE, row.names = F, quote = F)
#pca_umap<-read.table("/Users/vborda/Documents/Public_data/INS_data/PCA_UMAP/2024Jan_PCA_UMAP_INS_LDGH_1KGP_HGDP_nn50_MD0.5.txt", header=TRUE )

path<-"/Users/vborda/Documents/Public_data/INS_data/data+RefPops/"
demo<- read.table(paste0(path,"/FULL_1KGP_HGDP_PGP_IDs.txt"), head=FALSE, fill = TRUE,comment.char = "",sep = '\t')
key<- read.table(paste0(path,"UMAP/pop_colors.txt"), head=TRUE, fill = TRUE,comment.char = "")
colnames(demo)<-c("ID","POP","CodePop","City")
toplot<-merge(pca,demo, by.x = "ID", by.y = "ID")
#toplot<-merge(pca_umap,demo, by.x = "ID", by.y = "ID") 
toplot2<-merge(toplot,key, by.x = "POP", by.y = "POP")

##factor_variable <- factor(factor_variable, levels=c('this', 'that', 'those', ...))
toplot2$POP[toplot2$POP == "Ica"] <- "El Carmen"
toplot2$POP[toplot2$POP == "Ancash"] <- "Huaraz"
toplot2$POP[toplot2$POP == "Puno"] <- "Juliaca"
#toplot2$POP[toplot2$POP == "Ayacucho"] <- "Huamanga"
#toplot2$POP[toplot2$POP == "Moquegua"] <- "Mariscal Nieto"

#regions<-c("Tumbes","Lambayeque","Iquitos","Juliaca","Lima","El Carmen","Trujillo","Huaraz","Arequipa","Ayacucho",
#           "Moquegua","Cusco","Tacna")
regions<-c("PEL","Tumbes","Lambayeque","Iquitos","Juliaca","Lima","El Carmen","Trujillo","Huaraz","Arequipa","Ayacucho",
           "Moquegua","Cusco","Tacna")

references<-toplot2[!(toplot2$POP %in% regions),]
targets<-toplot2[(toplot2$POP %in% regions),]


inds_dataset<-as.data.frame(table(demo$POP))
inds_pca<-as.data.frame(table(toplot2$POP))

colnames(inds_dataset)<-c("POP","total dataset")
colnames(inds_pca)<-c("POP","total PCA")

df_merged <- merge(inds_dataset, inds_pca, by = "POP", all.x = TRUE)



library(ggpubr)
library(ggplot2)
library(gapminder)
library(gghalves)
library(ggdist)
library(ggforce)

#admixregions<-c("Iquitos","Tumbes","Lambayeque","Trujillo","Huaraz","Lima","El Carmen","Arequipa",
#                "Moquegua","Tacna","Ayacucho","Cusco","Juliaca")
admixregions<-c("PEL","Iquitos","Tumbes","Lambayeque","Trujillo","Huaraz","Lima","El Carmen","Arequipa",
                "Moquegua","Tacna","Ayacucho","Cusco","Juliaca")

#references<-c("GWD","MSL","YRI","LWK",
#              "CEU","TSI","IBS",
#              "CHB","CHS","JPT",
#              "Maya","Pima","Colombian",
#              "Tallanes","Moche","Quechuas","Chopccas","Qeros","Aimaras","Jacarus","Uros",
#              "Matses","Nahua","Shipibo","Lamas","Candoshi","Awajun","Matsiguenkas","Chachapoyas",
#              "Ashaninkas","Shimaa",
#              "PUR","MXL","CLM","PEL")

references<-c("GWD","MSL","YRI","LWK",
              "CEU","TSI","IBS",
              "CHB","CHS","JPT",
              "Maya","Pima","Colombian",
              "Tallanes","Moche","Quechuas","Chopccas","Qeros","Aimaras","Jacarus","Uros",
              "Matses","Nahua","Shipibo","Lamas","Candoshi","Awajun","Matsiguenkas","Chachapoyas",
              "Ashaninkas","Shimaa",
              "PUR","MXL","CLM")

KEY<-c(references,admixregions)

toplot2$POP<-factor(toplot2$POP,levels = KEY,ordered = TRUE)
scalesymbols <- toplot2$symbol
scalecolors <- toplot2$color


#jpeg(filename = "/Users/vborda/Documents/Public_data/INS_data/PCA_UMAP/20250117_PCA_1-2.jpg",
#     width = 22, height = 15, units = "cm", pointsize =12,  res = 600)

#ggplot(data = toplot2, aes(x = PC1, y = PC2,color=as.factor(POP),shape=as.factor(POP))) + 
#  geom_point(size=2,alpha = 0.5) +
#  scale_shape_manual(name="Populations",values = scalesymbols) +
#  scale_color_manual(name="Populations",values = scalecolors) +
#  theme(legend.text = element_text(size=8),
#        legend.title=element_text(size=8),
#        legend.position = c(0.2,0.65),legend.key.size = unit(0.45, 'cm'))

#dev.off()

###### OLD without PEL

#regions<-data.frame(
#  Deparment=c("Tumbes","Lambayeque","Iquitos","Trujillo","Ancash","Lima","Ica","Arequipa","Ayacucho","Moquegua","Tacna","Cusco","Puno"),
#  Province=c("Tumbes","Lambayeque","Iquitos","Trujillo","Huaraz","Lima","Chincha","Arequipa","Huamanga","Mariscal Nieto","Tacna",
#             "Cusco","Puno"),
#  Cities=c("Tumbes","Lambayeque","Iquitos","Trujillo","Huaraz","Lima","El Carmen","Arequipa","Ayacucho","Moquegua","Tacna",
#           "Cusco","Juliaca"),
#  tag=c(1:13),
#  symbol=c(7,21,8,22,9,23,10,24,12,25,13,4,14),
#  color=c("cornflowerblue","darkred","darkolivegreen1","darkorange","cyan","purple","gold1",
#          "cornflowerblue","darkred","darkolivegreen1","darkorange","cyan","purple")
#)

#REF<-data.frame(
#  Province=c("Tallanes","Moche","Chachapoyas","Awajun","Candoshi","Matses","Lamas","Ashaninkas","Matsiguenkas","Shimaa","Nahua","Colombian",
#             "Surui","Karitiana","Shipibo",
#             "Aimaras","Chopccas","Jacarus","Qeros","Quechuas","Uros",
#             "Pima","Maya","PUR","MXL","CLM","PEL",
#             "CEU","IBS","TSI",
#             "GWD","MSL","YRI","LWK",
#             "JPT","CHB","CHS"),
#  Area=c("Coast","Coast","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian",
#         "Amazonian","Amazonian","Amazonian",
#         "Andean","Andean","Andean","Andean","Andean","Andean",
#         "MesoAmerican","MesoAmerican","Admixed","Admixed","Admixed","Admixed",
#         "European","European","European",
#         "African","African","African","African",
#         "East Asian","East Asian","East Asian"),
#  symbol=c(4,4,rep(15, times = 13),
#           rep(16, times = 6),
#           11,11,6,6,6,6,
#           17,17,17,
#           18,18,18,18,
#           5,5,5)
#)
#REF$color<-"grey"


###### WITH PEL
regions<-data.frame(
  Deparment=c("PEL","Tumbes","Lambayeque","Iquitos","Trujillo","Ancash","Lima","Ica","Arequipa","Ayacucho","Moquegua","Tacna","Cusco","Puno"),
  Province=c("PEL","Tumbes","Lambayeque","Iquitos","Trujillo","Huaraz","Lima","Chincha","Arequipa","Huamanga","Mariscal Nieto","Tacna",
             "Cusco","Puno"),
  Cities=c("PEL","Tumbes","Lambayeque","Iquitos","Trujillo","Huaraz","Lima","El Carmen","Arequipa","Ayacucho","Moquegua","Tacna",
           "Cusco","Juliaca"),
  tag=c(1:14),
  symbol=c(15,7,21,8,22,9,23,10,24,12,25,13,4,14),
  color=c("tan","cornflowerblue","darkred","darkolivegreen1","darkorange","cyan","purple","gold1",
          "cornflowerblue","darkred","darkolivegreen1","darkorange","cyan","purple")
)


REF<-data.frame(
  Province=c("Tallanes","Moche","Chachapoyas","Awajun","Candoshi","Matses","Lamas","Ashaninkas","Matsiguenkas","Shimaa","Nahua","Colombian",
             "Surui","Karitiana","Shipibo",
             "Aimaras","Chopccas","Jacarus","Qeros","Quechuas","Uros",
             "Pima","Maya","PUR","MXL","CLM",
             "CEU","IBS","TSI",
             "GWD","MSL","YRI","LWK",
             "JPT","CHB","CHS"),
  Area=c("Coast","Coast","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian","Amazonian",
         "Amazonian","Amazonian","Amazonian",
         "Andean","Andean","Andean","Andean","Andean","Andean",
         "MesoAmerican","MesoAmerican","Admixed","Admixed","Admixed",
         "European","European","European",
         "African","African","African","African",
         "East Asian","East Asian","East Asian"),
  symbol=c(4,4,rep(15, times = 13),
           rep(16, times = 6),
           11,11,6,6,6,
           17,17,17,
           18,18,18,18,
           5,5,5)
)
REF$color<-"grey"




references<-toplot2[!(toplot2$POP %in% regions$Province),]
references$alpha<-0.1
#references$symbol<-19

#targets<-toplot2[(toplot2$POP %in% regions$Province),]
targets<-toplot2[(toplot2$POP %in% regions$Cities),]
targets$alpha<-1
#targets<-merge(targets[,c(1:12)],regions,by.x="POP",by.y="Province")
targets<-merge(targets[,c(1:12)],regions,by.x="POP",by.y="Cities")
references<-merge(references[,c(1:12)],REF,by.x="POP",by.y="Province")

library(tidyverse)
library(palmerpenguins)

targets$symbol<-as.factor(targets$symbol)
#targets$POP<-factor(targets$POP,
#                       levels=c("Tumbes","Lambayeque","Iquitos","Trujillo","Huaraz","Lima","El Carmen","Arequipa","Ayacucho",
#                                "Moquegua","Tacna","Cusco","Juliaca"))
targets$POP<-factor(targets$POP,
                    levels=c("PEL","Tumbes","Lambayeque","Iquitos","Trujillo","Huaraz","Lima","El Carmen","Arequipa","Ayacucho",
                             "Moquegua","Tacna","Cusco","Juliaca"))


p1<-ggplot() +              
  geom_point(data=references,aes(PC1,PC2),shape=references$symbol,alpha=0.5,color=references$color,size=3) +
  geom_point(data=targets,aes(PC1,PC2,shape=POP,color=POP),alpha=1,size=2)+
  scale_shape_manual(name = "Urban Population", 
                     values = c("PEL"=15,"Tumbes" = 7, "Lambayeque" = 21, "Iquitos" = 8, "Trujillo" = 22, "Huaraz" = 9,"Lima" = 23, "El Carmen" = 10, 
                    #values = c("Tumbes" = 7, "Lambayeque" = 21, "Iquitos" = 8, "Trujillo" = 22, "Huaraz" = 9,"Lima" = 23, "El Carmen" = 10, 
                                "Arequipa" = 24,"Ayacucho" = 12, "Moquegua" = 25,"Tacna" = 13, "Cusco" = 4, 
                                "Juliaca" = 14)) +
  scale_color_manual(name = "Urban Population", 
                     values = c("PEL"="tan","Tumbes" = "cornflowerblue","Lambayeque" = "darkred", "Iquitos" = "darkolivegreen1", 
                    #values = c("Tumbes" = "cornflowerblue","Lambayeque" = "darkred", "Iquitos" = "darkolivegreen1", 
                                "Trujillo" = "darkorange", "Huaraz" = "cyan", "Lima" = "purple", "El Carmen" = "gold1", 
                                "Arequipa" = "cornflowerblue", "Ayacucho" = "darkred", "Moquegua" = "darkolivegreen1", 
                                "Tacna" = "darkorange", "Cusco" = "cyan", "Juliaca" = "purple")) +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent'),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9, lineheight = 0.5),  # Reduces space between words in the text
    legend.spacing.y = unit(0.18, "cm") , # Reduces vertical space between items in the legend
    legend.key.height = unit(0.5, "cm"),   # Height of legend boxes
    legend.key.width = unit(0.5, "cm")
  )


p1_noL<-ggplot() +              
  geom_point(data=references,aes(PC1,PC2),shape=references$symbol,alpha=0.5,color=references$color,size=3) +
  geom_point(data=targets,aes(PC1,PC2,shape=POP,color=POP),alpha=1,size=2)+
  scale_shape_manual(name = "Urban Population", 
                     values = c("PEL"=15,"Tumbes" = 7, "Lambayeque" = 21, "Iquitos" = 8, "Trujillo" = 22, "Huaraz" = 9,"Lima" = 23, "El Carmen" = 10, 
                                #values = c("Tumbes" = 7, "Lambayeque" = 21, "Iquitos" = 8, "Trujillo" = 22, "Huaraz" = 9,"Lima" = 23, "El Carmen" = 10, 
                                "Arequipa" = 24,"Ayacucho" = 12, "Moquegua" = 25,"Tacna" = 13, "Cusco" = 4, 
                                "Juliaca" = 14)) +
  scale_color_manual(name = "Urban Population", 
                     values = c("PEL"="tan","Tumbes" = "cornflowerblue","Lambayeque" = "darkred", "Iquitos" = "darkolivegreen1", 
                                #values = c("Tumbes" = "cornflowerblue","Lambayeque" = "darkred", "Iquitos" = "darkolivegreen1", 
                                "Trujillo" = "darkorange", "Huaraz" = "cyan", "Lima" = "purple", "El Carmen" = "gold1", 
                                "Arequipa" = "cornflowerblue", "Ayacucho" = "darkred", "Moquegua" = "darkolivegreen1", 
                                "Tacna" = "darkorange", "Cusco" = "cyan", "Juliaca" = "purple")) +
  theme(plot.margin = margin(0.4,0.4,0.4,0.4, "cm"),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="none"
  )


p2_noL<-ggplot() +              
  geom_point(data=references,aes(PC3,PC4),shape=references$symbol,alpha=0.5,color=references$color,size=3) +
  geom_point(data=targets,aes(PC3,PC4,shape=POP,color=POP),alpha=1,size=2)+
  scale_shape_manual(name = "Urban Population", 
                     values = c("PEL"=15,"Tumbes" = 7, "Lambayeque" = 21, "Iquitos" = 8, "Trujillo" = 22, "Huaraz" = 9,"Lima" = 23, "El Carmen" = 10, 
                                #values = c("Tumbes" = 7, "Lambayeque" = 21, "Iquitos" = 8, "Trujillo" = 22, "Huaraz" = 9,"Lima" = 23, "El Carmen" = 10, 
                                "Arequipa" = 24,"Ayacucho" = 12, "Moquegua" = 25,"Tacna" = 13, "Cusco" = 4, 
                                "Juliaca" = 14)) +
  scale_color_manual(name = "Urban Population", 
                     values = c("PEL"="tan","Tumbes" = "cornflowerblue","Lambayeque" = "darkred", "Iquitos" = "darkolivegreen1", 
                                #values = c("Tumbes" = "cornflowerblue","Lambayeque" = "darkred", "Iquitos" = "darkolivegreen1", 
                                "Trujillo" = "darkorange", "Huaraz" = "cyan", "Lima" = "purple", "El Carmen" = "gold1", 
                                "Arequipa" = "cornflowerblue", "Ayacucho" = "darkred", "Moquegua" = "darkolivegreen1", 
                                "Tacna" = "darkorange", "Cusco" = "cyan", "Juliaca" = "purple")) +
  theme(plot.margin = margin(0.4,0.4,0.4,0.4, "cm"),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="none"
  )

p3_noL<-ggplot() +              
  geom_point(data=references,aes(PC5,PC6),shape=references$symbol,alpha=0.5,color=references$color,size=3) +
  geom_point(data=targets,aes(PC5,PC6,shape=POP,color=POP),alpha=1,size=2)+
  scale_shape_manual(name = "Urban Population", 
                     values = c("PEL"=15,"Tumbes" = 7, "Lambayeque" = 21, "Iquitos" = 8, "Trujillo" = 22, "Huaraz" = 9,"Lima" = 23, "El Carmen" = 10, 
                                #values = c("Tumbes" = 7, "Lambayeque" = 21, "Iquitos" = 8, "Trujillo" = 22, "Huaraz" = 9,"Lima" = 23, "El Carmen" = 10, 
                                "Arequipa" = 24,"Ayacucho" = 12, "Moquegua" = 25,"Tacna" = 13, "Cusco" = 4, 
                                "Juliaca" = 14)) +
  scale_color_manual(name = "Urban Population", 
                     values = c("PEL"="tan","Tumbes" = "cornflowerblue","Lambayeque" = "darkred", "Iquitos" = "darkolivegreen1", 
                                #values = c("Tumbes" = "cornflowerblue","Lambayeque" = "darkred", "Iquitos" = "darkolivegreen1", 
                                "Trujillo" = "darkorange", "Huaraz" = "cyan", "Lima" = "purple", "El Carmen" = "gold1", 
                                "Arequipa" = "cornflowerblue", "Ayacucho" = "darkred", "Moquegua" = "darkolivegreen1", 
                                "Tacna" = "darkorange", "Cusco" = "cyan", "Juliaca" = "purple")) +
  theme(plot.margin = margin(0.4,0.4,0.4,0.4, "cm"),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="none"
  )


p4_noL<-ggplot() +              
  geom_point(data=references,aes(PC7,PC8),shape=references$symbol,alpha=0.5,color=references$color,size=3) +
  geom_point(data=targets,aes(PC7,PC8,shape=POP,color=POP),alpha=1,size=2)+
  scale_shape_manual(name = "Urban Population", 
                     values = c("PEL"=15,"Tumbes" = 7, "Lambayeque" = 21, "Iquitos" = 8, "Trujillo" = 22, "Huaraz" = 9,"Lima" = 23, "El Carmen" = 10, 
                                #values = c("Tumbes" = 7, "Lambayeque" = 21, "Iquitos" = 8, "Trujillo" = 22, "Huaraz" = 9,"Lima" = 23, "El Carmen" = 10, 
                                "Arequipa" = 24,"Ayacucho" = 12, "Moquegua" = 25,"Tacna" = 13, "Cusco" = 4, 
                                "Juliaca" = 14)) +
  scale_color_manual(name = "Urban Population", 
                     values = c("PEL"="tan","Tumbes" = "cornflowerblue","Lambayeque" = "darkred", "Iquitos" = "darkolivegreen1", 
                                #values = c("Tumbes" = "cornflowerblue","Lambayeque" = "darkred", "Iquitos" = "darkolivegreen1", 
                                "Trujillo" = "darkorange", "Huaraz" = "cyan", "Lima" = "purple", "El Carmen" = "gold1", 
                                "Arequipa" = "cornflowerblue", "Ayacucho" = "darkred", "Moquegua" = "darkolivegreen1", 
                                "Tacna" = "darkorange", "Cusco" = "cyan", "Juliaca" = "purple")) +
  theme(plot.margin = margin(0.4,0.4,0.4,0.4, "cm"),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="none"
  )

p5_noL<-ggplot() +              
  geom_point(data=references,aes(PC9,PC10),shape=references$symbol,alpha=0.5,color=references$color,size=3) +
  geom_point(data=targets,aes(PC9,PC10,shape=POP,color=POP),alpha=1,size=2)+
  scale_shape_manual(name = "Urban Population", 
                     values = c("PEL"=15,"Tumbes" = 7, "Lambayeque" = 21, "Iquitos" = 8, "Trujillo" = 22, "Huaraz" = 9,"Lima" = 23, "El Carmen" = 10, 
                                #values = c("Tumbes" = 7, "Lambayeque" = 21, "Iquitos" = 8, "Trujillo" = 22, "Huaraz" = 9,"Lima" = 23, "El Carmen" = 10, 
                                "Arequipa" = 24,"Ayacucho" = 12, "Moquegua" = 25,"Tacna" = 13, "Cusco" = 4, 
                                "Juliaca" = 14)) +
  scale_color_manual(name = "Urban Population", 
                     values = c("PEL"="tan","Tumbes" = "cornflowerblue","Lambayeque" = "darkred", "Iquitos" = "darkolivegreen1", 
                                #values = c("Tumbes" = "cornflowerblue","Lambayeque" = "darkred", "Iquitos" = "darkolivegreen1", 
                                "Trujillo" = "darkorange", "Huaraz" = "cyan", "Lima" = "purple", "El Carmen" = "gold1", 
                                "Arequipa" = "cornflowerblue", "Ayacucho" = "darkred", "Moquegua" = "darkolivegreen1", 
                                "Tacna" = "darkorange", "Cusco" = "cyan", "Juliaca" = "purple")) +
  theme(plot.margin = margin(0.4,0.4,0.4,0.4, "cm"),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="none"
  )

#### arrows

arrows_p1 <- 
  tibble(
    x1 = c(-0.04,-0.012,0.02,0.011,-0.024,-0.024),
    x2 = c(-0.04,-0.002,0.0076,0.022,-0.021,-0.007),
    y1 = c(-0.007,0.04,0.009,-0.025,0.017,0.017),
    y2 = c(-0.015,0.04,0.0047,-0.025,0.009,0.0158)
  )

P1<-p1_noL +  
  annotate("text", x = -0.04, y = -0.005, label = "AFRICAN", colour = "blue", fontface =2,size=2) +
  annotate("text", x = -0.018, y = 0.04, label = "EUROPEAN", colour = "darkred", fontface =2,size=2)+
  annotate("text", x = 0.02, y = 0.01, label = "EAST ASIAN", colour = "gold4", fontface =2,size=2)+
  annotate("text", x = 0, y = -0.025, label = "INDIGENOUS AMERICAN", colour = "darkgreen", fontface =2,size=2)+
  annotate("text", x = -0.025, y = 0.02, label = "1KGP\nLATIN AMERICAN", colour = "grey15", fontface =2,size=2)+
  geom_curve(
    data = arrows_p1, aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.08, "inch"),type="closed"), size = 0.5,
    color = "gray20", curvature = 0)

arrows_p2 <- 
  tibble(#LATAM,ANDEAN, AMAZONIAN
    x1 = c(0.025,0.002,-0.001),
    x2 = c(0.016,-0.008,-0.009),
    y1 = c(0.003,-0.03,0.1),
    y2 = c(-0.007,-0.035,0.1)
  )

P2<-p2_noL +  
  annotate("text", x = 0.048, y = 0.013, label = "EAST ASIAN", colour = "gold4", fontface =2,size=2)+
  annotate("text", x = 0.008, y = -0.028, label = "INDIGENOUS AMERICAN\n(ANDEAN)", colour = "tan", fontface =2,size=2,vjust=0.5)+
  annotate("text", x = 0.008, y = 0.1, label = "INDIGENOUS AMERICAN\n(AMAZONIAN)", colour = "darkgreen", fontface =2,size=2,vjust=0.5)+
  annotate("text", x = 0.025, y = 0.01, label = "1KGP\nLATIN AMERICAN", colour = "grey15", fontface =2,size=2)+
  geom_curve(
    data = arrows_p2, aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.08, "inch"),type="closed"), size = 0.5,
    color = "gray20", curvature = 0)


arrows_p3 <- 
  tibble(
    x1 = c(0.05,-0.08,-0.05,-0.05),
    x2 = c(0.05,-0.08,0.0076,-0.0052),
    y1 = c(0.013,0.028,0.11,-0.047),
    y2 = c(0,0.01,0.11,-0.065)
  )

P3<-p3_noL +  
  annotate("text", x = 0.053, y = 0.02, label = "AFRICAN\n(GWD)", colour = "blue", fontface =2,size=2,vjust=0.5) +
  annotate("text", x = -0.08, y = 0.035, label = "AFRICAN\n(LWK)", colour = "blue", fontface =2,size=2,vjust=0.5) +
  annotate("text", x = -0.05, y = 0.12, label = "INDIGENOUS AMERICAN\n(SOUTH AMAZONIAN-Shimaa)", colour = "darkgreen", fontface =2,size=2,vjust=0.5)+
  annotate("text", x = -0.05, y = -0.04, label = "INDIGENOUS AMERICAN\n(NORTH AMAZONIAN-Awajun)", colour = "darkgreen", fontface =2,size=2,vjust=0.5)+
  geom_curve(
    data = arrows_p3, aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.08, "inch"),type="closed"), size = 0.5,
    color = "gray20", curvature = 0)


P4<-p4_noL +  
  annotate("text", x = -0.04, y = -0.05, label = "EAST ASIAN\n(CHS)", colour = "gold4", fontface =2,size=2,vjust=0.5) +
  annotate("text", x = 0.05, y = 0.05, label = "EAST ASIAN\n(JPT)", colour = "gold4", fontface =2,size=2,vjust=0.5) +
  annotate("text", x = -0.005, y = 0.105, label = "INDIGENOUS AMERICAN\n(SOUTH AMAZONIAN-Shimaa)", colour = "darkgreen", fontface =2,size=2,vjust=0.5)+
  annotate("text", x = 0.045, y = -0.1, label = "INDIGENOUS AMERICAN\n(Pima)", colour = "darkgreen", fontface =2,size=2,vjust=0.5)

arrows_p5 <- 
  tibble(
    x1 = c(-0.07,0.065,-0.05,-0.05),
    x2 = c(-0.07,0.065,-0.01,0),
    y1 = c(-0.018,0.028,-0.22,0.085),
    y2 = c(0,0.01,-0.22,0.085)
  )

P5<-p5_noL +  
  annotate("text", x = -0.07, y = -0.03, label = "AFRICAN\n(GWD)", colour = "blue", fontface =2,size=2,vjust=0.5) +
  annotate("text", x = 0.065, y = 0.041, label = "AFRICAN\n(YRI)", colour = "blue", fontface =2,size=2,vjust=0.5) +
  annotate("text", x = -0.05, y = -0.2, label = "INDIGENOUS AMERICAN\n(Pima)", colour = "darkgreen", fontface =2,size=2,vjust=0.5)+
  annotate("text", x = -0.05, y = 0.1, label = "INDIGENOUS AMERICAN\n(SOUTH AMAZONIAN-Shimaa)", colour = "darkgreen", fontface =2,size=2,vjust=0.5)+
  geom_curve(
    data = arrows_p5, aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.08, "inch"),type="closed"), size = 0.5,
    color = "gray20", curvature = 0)



library("ggplot2") 
library("grid") 
library("gridExtra") 
library("cowplot") 

legend <- get_legend(p1) 
grid.newpage()
grid.draw(legend)  

library(ggpubr)
jpeg(filename = "/Users/vborda/Documents/Public_data/INS_data/Manuscript/R1_Figures/Fig_S2_20250808_PCA_1-10.jpg",
#jpeg(filename = "/Users/vborda/Documents/Public_data/INS_data/PCA_UMAP/Fig_S2_20250117_PCA_1-10.jpg",
     width = 21, height = 25, units = "cm", pointsize =12,  res = 600)


ggarrange(P1,P2,P3,P4,P5,legend, ncol=2,nrow=3, labels = c("A","B","C","D","E"), vjust = 2.2,
                        font.label=list(color="black",size=11))

dev.off()



