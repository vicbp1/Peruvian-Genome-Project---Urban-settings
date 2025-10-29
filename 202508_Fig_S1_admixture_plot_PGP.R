
#### ADMIXTURE results
## Extracting K results

## for K in $(seq 3 10); do grep ^Loglike log_K${K}_run* | sort -gk2 | tail -1 \
## | cut -d":" -f1 | while read file; do grep CV $file ; done ; done | sed -e "s/CV error (//g" -e "s/): /\t/g" > 


## paste POP_ID_input2admixture outputs/K10_run_7.Q > proportions/admixture_proportions_K10.txt


autopropnamefile<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_K4.txt"
xpropnamefile<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_chrX_K4.txt"

autopropnamefile5<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_K5.txt"
xpropnamefile5<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_chrX_K5.txt"

autopropnamefile6<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_K6.txt"
autopropnamefile7<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_K7.txt"
autopropnamefile8<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_K8.txt"
autopropnamefile9<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_K9.txt"
autopropnamefile10<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_K10.txt"
labels_for_plot<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/labels_for_admixture"


autoprop4<-read.table(autopropnamefile,header=FALSE,sep="")
xprop4<-read.table(xpropnamefile,header=FALSE,sep="\t")

autoprop5<-read.table(autopropnamefile5,header=FALSE,sep="")
xprop5<-read.table(xpropnamefile5,header=FALSE,sep="\t")

autoprop6<-read.table(autopropnamefile6,header=FALSE,sep="")
autoprop7<-read.table(autopropnamefile7,header=FALSE,sep="")
autoprop8<-read.table(autopropnamefile8,header=FALSE,sep="")
autoprop9<-read.table(autopropnamefile9,header=FALSE,sep="")
autoprop10<-read.table(autopropnamefile10,header=FALSE,sep="")
labels_plot<-read.table(labels_for_plot,header=FALSE,sep="\t")

colnames(labels_plot)<-c("Population","id","AFR","EUR","EAS","NAT","Continental Groups")
TwoLabels<-labels_plot[,c(1,7)]
##factor_variable <- factor(factor_variable, levels=c('this', 'that', 'those', ...))
TwoLabels$Population[TwoLabels$Population == "Afro_des"] <- "El Carmen"
TwoLabels$Population[TwoLabels$Population == "Ancash"] <- "Huaraz"
TwoLabels$Population[TwoLabels$Population == "Puno"] <- "Juliaca"
#TwoLabels$Population[TwoLabels$Population == "Ayacucho"] <- "Huamanga"
#TwoLabels$Population[TwoLabels$Population == "Moquegua"] <- "Mariscal Nieto"

colnames(autoprop4)<-c("Population","id","AFR","EUR","EAS","NAT")
colnames(autoprop5)<-c("Population","id","EAS","EUR","AFR","ANDEAN","AMAZONIAN")
#colnames(xprop4)<-c("Population","id","AFR","EUR","EAS","NAT")
#colnames(xprop5)<-c("Population","id","AFR","EUR","EAS","NAT")
#colnames(autoprop6)<-c("Population","id","Amazonian","African","South Andean","North Andean","East Asian","European")
#colnames(xprop6)<-c("Population","id","African","EUR1","EAS","Amazonian","EUR2",)


library(xadmix)  
library(ggplot2)  
library(gridExtra)  
library(label.switching)  
library(tidyr)  
library(remotes)

#remotes::install_github('royfrancis/pophelper')
library(pophelper)

## Creating fake 10 for i in $(seq 1 2162); do echo -e "0\t0\t0\t0\t0\t0\t0\t0\t0\t0" ; done > K10.txt

afiles <- list.files(path="/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/topophelper",
                     full.names=TRUE)
q1 <- readQ(files=afiles)
# create a qlist

#q1 <- list("K4"=autoprop4[,2:6],"K5"=autoprop5[,2:7])
#str(q1)



autoprop5$Population[autoprop5$Population == "Afro_des"] <- "El Carmen"
autoprop5$Population[autoprop5$Population == "Ancash"] <- "Huaraz"
autoprop5$Population[autoprop5$Population == "Puno"] <- "Juliaca"
autoprop5$Population[autoprop5$Population == "Shipibo_INS"] <- "Shipibo"
#autoprop4$Population[autoprop4$Population == "Ayacucho"] <- "Huamanga"
#autoprop4$Population[autoprop4$Population == "Moquegua"] <- "Mariscal Nieto"

admixregions<-c("Iquitos","Tumbes","Lambayeque","Trujillo","Huaraz","Lima","El Carmen","Arequipa","Moquegua","Tacna","Ayacucho","Cusco","Juliaca")
references<-c("GWD","MSL","YRI","LWK",
              "CEU","TSI","IBS",
              "CHB","CHS","JPT",
              "Maya","Pima","Colombian",
              "Tallanes","Moche","Quechuas","Chopccas","Qeros","Aimaras","Jacarus","Uros",
              "Matses","Nahua","Shipibo","Lamas","Candoshi","Awajun","Matsiguenkas","Chachapoyas",
              "Ashaninkas","Shimaa",
              "PUR","MXL","CLM","PEL")
allpops=c(admixregions,references)

colnames(autoprop5)[colnames(autoprop5) == "Population"] ="Pop"


#plotQ(qlist=q1[c(1:4)], grplab=autoprop4[,1,drop=FALSE],exportpath="/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/",
#      subsetgrp=admixregions,imgoutput="join",basesize=2.5, ## basezise # all names outside the plot
#      grplabspacer=0.05,grplabsize=0.8,grplabpos=1,grplabjust=0.4,pointsize=6,linesize=7,linealpha=0.2,  ## pointsize: linea blanca
#      pointcol="white",grplabangle=-90,linepos=0.8,grplabheight=0.5,sharedindlab=TRUE,
#      divcol="black",divtype=1,divsize=0.3,
#      showlegend=FALSE,legendkeysize=5,legendtextsize=3,legendlab=c("African 1","European","Andean 1","East Asian","Amazonian","Andean 2"),
#      panelspacer = 0.1,splab = c("K=3","K=4","K=5","K=6"), splabsize=2,showyaxis=FALSE,
#      clustercol=c("darkblue","firebrick2","burlywood2","yellow1","darkgreen","lightgreen","dodgerblue3"),ordergrp=TRUE, dpi=600)


plot_Peru<-plotQ(qlist=q1[c(1:7)], grplab=autoprop5[,1,drop=FALSE],returnplot = TRUE,exportplot=F ,#exportpath="/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/",
      subsetgrp=admixregions,imgoutput="join",basesize=2.5, showtitle = TRUE,titlelab = "B", titlesize = 12, titleface = "bold", ## basezise # all names outside the plot
      grplabspacer=0.05,grplabsize=1.8,grplabpos=0.25,grplabjust=1,pointsize=0,linesize=0,linealpha=0.2,  ## pointsize: linea blanca
      pointcol="black",grplabangle=-90,linepos=0.8,grplabheight=20,sharedindlab=TRUE,pointbgcol = "white",# grplabbgcol="white",
      divcol="black",divtype=1,divsize=0.3,indlabheight = 1,
      showlegend=FALSE,
      panelspacer = 0.1,splab = c("K=3","K=4","K=5","K=6","K=7","K=8","K=9"), splabsize=7,showyaxis=FALSE,
      clustercol=c("darkblue","firebrick2","burlywood2","yellow1","darkgreen","lightgreen","dodgerblue3","purple","gold4"),ordergrp=TRUE, dpi=600)


plot_Peru2<-plotQ(qlist=q1[c(1:8)], grplab=autoprop5[,1,drop=FALSE],returnplot = TRUE,exportplot=F ,#exportpath="/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/",
                 subsetgrp=admixregions,imgoutput="join",basesize=2.5, showtitle = TRUE,titlelab = "B", titlesize = 12, titleface = "bold", ## basezise # all names outside the plot
                 grplabspacer=0.05,grplabsize=1.8,grplabpos=0.25,grplabjust=1,pointsize=0,linesize=0,linealpha=0.2,  ## pointsize: linea blanca
                 pointcol="black",grplabangle=-90,linepos=0.8,grplabheight=20,sharedindlab=TRUE,pointbgcol = "white",# grplabbgcol="white",
                 divcol="black",divtype=1,divsize=0.3,indlabheight = 1,
                 showlegend=FALSE,
                 panelspacer = 0.1,splab = c("K=3","K=4","K=5","K=6","K=7","K=8","K=9","K=10"), splabsize=7,showyaxis=FALSE,
                 clustercol=c("darkblue","firebrick2","burlywood2","yellow1","darkgreen","lightgreen","dodgerblue3","purple","gold4","turquoise2"),ordergrp=TRUE, dpi=600)


### ALL

#plot_ALL<-plotQ(qlist=q1[c(2:8)], grplab=autoprop4[,1,drop=FALSE],returnplot = TRUE,exportplot=F, #exportpath="/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/",
#      subsetgrp=allpops,imgoutput="join",basesize=2.5, ## basezise # all names outside the plot
#      grplabspacer=0.05,grplabsize=0.8,grplabpos=1,grplabjust=0.3,pointsize=8,linesize=7,linealpha=0.2,  ## pointsize: linea blanca
#      pointcol="white",grplabangle=-90,linepos=0.8,grplabheight=2,sharedindlab=TRUE,
#      divcol="black",divtype=1,divsize=0.3,
#      showlegend=FALSE,
#      panelspacer = 0.1,splab = c("K=3","K=4","K=5","K=6","K=7","K=8","K=9"), splabsize=20,showyaxis=FALSE,
#      clustercol=c("darkblue","firebrick2","burlywood2","yellow1","darkgreen","lightgreen","dodgerblue3","purple","gold4"),ordergrp=TRUE, dpi=600)




plot_REF<-plotQ(qlist=q1[c(1:7)], grplab=autoprop4[,1,drop=FALSE],,returnplot = TRUE,exportplot=F, #exportpath="/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/",
                subsetgrp=references,imgoutput="join",basesize=1, showtitle = TRUE,titlelab = "A" , titlesize = 12, titleface = "bold",## basezise # all names outside the plot
                grplabheight=1,grplabspacer=0.05,grplabsize=1.8,grplabpos=0.25,grplabjust=1,pointsize=0,linesize=0,linealpha=0.2,  ## pointsize: linea blanca ,grplabsize=3, line size =shadow
                pointcol="black",grplabangle=-90,linepos=1,sharedindlab=TRUE,pointbgcol = "white",
                divcol="black",divtype=1,divsize=0.3,indlabheight = 1,
                showlegend=FALSE,
                panelspacer = 0.1,splab = c("K=3","K=4","K=5","K=6","K=7","K=8","K=9"), splabsize=7,showyaxis=FALSE,
                clustercol=c("darkblue","firebrick2","burlywood2","yellow1","darkgreen","lightgreen","dodgerblue3","purple","gold4"),ordergrp=TRUE, dpi=600)


plot_REF2<-plotQ(qlist=q1[c(1:8)], grplab=autoprop5[,1,drop=FALSE],,returnplot = TRUE,exportplot=F, #exportpath="/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/",
                subsetgrp=references,imgoutput="join",basesize=1, showtitle = TRUE,titlelab = "A" , titlesize = 12, titleface = "bold",## basezise # all names outside the plot
                grplabheight=1,grplabspacer=0.05,grplabsize=1.8,grplabpos=0.25,grplabjust=1,pointsize=0,linesize=0,linealpha=0.2,  ## pointsize: linea blanca ,grplabsize=3, line size =shadow
                pointcol="black",grplabangle=-90,linepos=1,sharedindlab=TRUE,pointbgcol = "white",
                divcol="black",divtype=1,divsize=0.3,indlabheight = 1,
                showlegend=FALSE,
                panelspacer = 0.1,splab = c("K=3","K=4","K=5","K=6","K=7","K=8","K=9","K=10"), splabsize=7,showyaxis=FALSE,
                clustercol=c("darkblue","firebrick2","burlywood2","yellow1","darkgreen","lightgreen","dodgerblue3","purple","gold4","turquoise2"),ordergrp=TRUE, dpi=600)



library(ggpubr)


          #labels = c("PC1 vs PC2", "PC3 vs PC4", "PC5 vs PC6", "PC7 vs PC8", "PC9 vs PC10"),
          #ncol = 2, nrow = 3)




## K10


#plot_Peru<-plotQ(qlist=q1[1], grplab=autoprop4[,1,drop=FALSE],returnplot = TRUE,exportpath="/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/",
#                 subsetgrp=admixregions,basesize=2.5, ## basezise # all names outside the plot
#                 grplabspacer=0.05,grplabsize=0.8,grplabpos=1,grplabjust=0.4,pointsize=1,linesize=7,linealpha=0.2,  ## pointsize: linea blanca
#                 pointcol="white",grplabangle=90,linepos=0.8,grplabheight=0.5,sharedindlab=TRUE,
#                 divcol="black",divtype=1,divsize=0.3,indlabheight = 1,
#                 showlegend=FALSE,
#                 panelspacer = 0.1,splab = c("K10"), splabsize=2,showyaxis=FALSE,
#                 clustercol=c("darkblue","firebrick2","burlywood2","yellow1","darkgreen","lightgreen","dodgerblue3","purple","cadetblue","blue"),ordergrp=TRUE, dpi=600)




CV=read.table("/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/K_CV_values.txt",header=FALSE,sep="\t")
colnames(CV)=c("K value","Cross-validation value")

CV$`K value`=factor(CV$`K value`,levels=CV$`K value`,ordered=TRUE)

CV_plot<-ggplot(CV, aes(x = `K value`, y = `Cross-validation value`)) +  geom_line(group=1,color="firebrick") +
  geom_point(color = "firebrick") +
  theme(plot.margin = margin(2,0.1,2,0.1, "cm"),
        axis.text=element_text(size=6),
        axis.title=element_text(size=7,face = "bold"))


library(ggpubr)

bottom_row <- ggarrange(plot_Peru2$plot[[1]],CV_plot, ncol=2,nrow=1,widths = c(1,0.5), labels = c(" ","C"), vjust = 2.2,
                        font.label=list(color="black",size=11))

#garrange(p1,bottom_row,ncol=1,nrow=2)

labels = c("A", "D", "B", "F", "C", "G", "", "")

#ggsave("/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/20250107_ADMIXTURE_PGP.jpg",
#       grid.arrange(plot_REF$plot[[1]],,ncol=1,heights=c(1,1,0.5))
#       , width =18, height = 21, units = "cm") # Figure 4A


ggsave("/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/Fig_S1_20250824_ADMIXTURE_PGP.jpg",
       ggarrange(plot_REF2$plot[[1]],bottom_row,ncol=1,heights=c(1,1)),
       width =18, height = 22, units = "cm") # Figure 4A




#plotQ(qlist=q1[1:2],imgoutput="join",exportpath="/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/")



k3plot <-
  ggplot(kdf3, aes(factor(sampleID), prob, fill = factor(popGroup))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(loc), scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=3", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE)