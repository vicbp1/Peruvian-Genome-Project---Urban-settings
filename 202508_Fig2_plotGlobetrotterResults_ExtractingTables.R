

### SCRIPT TO PLOT GLOBETROTTER RESULTS FOR PAPER ###

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################

main_dir <-"/Users/vborda/Documents/Public_data/INS_data/Chromopainter_runs/globetrotter/globetrotter_Nov2024/Results/"
GT_dir <-"/Users/vborda/Documents/general_pipelines/plotting_scripts/globetrotter_scripts/"
plots<-"/Users/vborda/Documents/Public_data/INS_data/Chromopainter_runs/globetrotter/globetrotter_2025/"

source(paste0(GT_dir,"getLines.R"))
source(paste0(GT_dir,"makeDate_Wang2023.R"))
source(paste0(GT_dir,"getGlobetrotter.R"))
#source(paste0(GT_dir,"getGlobetrotterDates_onefile.R"))
#source(paste0(GT_dir,"getFastGlobetrotterDates.R"))
source(paste0(GT_dir,"addGTresults.R"))

##IMPORTANTE EN tabla TO PLOT CAMBIAS LOS NA de R por NA reales
##########################################################
## DEFINE DATAFILES

res_dir <- main_dir              	## Results directory
popkey_file <- "/Users/vborda/Documents/Public_data/INS_data/Chromopainter_runs/globetrotter/PGP_PopulationKey_500K_August2025.txt"
popkey <- read.table(popkey_file,header=T)
popkey$Ethnic_Group <- toupper(popkey$Population)                     	## Add a column
popplot <- as.character(popkey$Population)                           	## creating a character variable of population names

leginfo <- read.table(popkey_file, header = T, comment.char = "")
#ancreg_list <- c("MSL","YRI","Ghana","Herero","Mbukushu","LWK","Uganda" )
ancreg_list <- as.character(popplot[14:length(popplot)])

############ IMPORTANTE PARA EL ORDEN DE LAS POBLACIONES A PLOTAR
popplot2 <- popplot[1:14]                               ### 9 indicates the number of populations to plot
popplot <- popplot[15:length(popplot)]
popplotorder <- popplot
popplot <- factor(popplot,levels=popplotorder)
popplotorder2 <- popplot2
popplot2 <- factor(popplot2,levels=popplotorder2)

## LOAD POPKEY FILE ##

### DEFINE SOME PLOTTING VARIABLES ##
dateLabelCex=1
#datelines=c(1600,1650,1700,1750,1800,1850,1900,1950,2000)
#yAxisLim=c(1600,2000)

#datelines=c(1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80)
datelines=c(1,5,10,15,20)
yAxisLim=c(1,20)

###############################################################################


leginfo2order <- leginfo[match(popplot,leginfo$Population),]

pcolshex <-leginfo2order$Colour
## PULL IN DATA

## PULL IN THE RESULTS OF THE MIXTURE MODEL (RUN ELSEWHERE)

inputname<-"globetrotter_results_Admixed_PGP_20250812"

final.res2plot <- read.table(paste0(main_dir,inputname,"_toplot.txt"),header=T,as.is=T)
#mixmat <- read.table("/Users/vborda/Documents/Public_data/INS_data/Chromopainter_runs/contributions/Admixed_1.9M_sourcefind_2024.txt",header=T,row.names=1,as.is=T)
admixturesources2 <- read.table(paste0(main_dir,inputname,"_AdmixtureSources2.txt"),header=T,as.is=T)
admixturesources2 <- admixturesources2[,c(14:ncol(admixturesources2))]
dateboots <- read.table(paste0(main_dir,inputname,"_OneDateBootstraps.txt"),header=T,row.names=1,as.is=T)
date2boots <- read.table(paste0(main_dir,inputname,"_TwoDateBootstraps.txt"),header=T,row.names=1,as.is=T)


## Extracting data table admixture sources:

admixturesources_subset <- admixturesources2[grep("null\\.(1|2)$", rownames(admixturesources2)), ]
rownames(admixturesources_subset) <- sub("null", "source", rownames(admixturesources_subset))
admixturesources_subset <- admixturesources_subset[, colSums(admixturesources_subset) != 0]

library(tidyverse)

admixturesources_subset <- admixturesources_subset %>%
  rownames_to_column("rowname") %>%                                  # move row names into a column
  separate(rowname, into = c("city", "rest"), sep = "\\.", extra = "merge") %>%  # split at first dot
  separate(rest, into = c("source", "id"), sep = "\\.", fill = "right")          # split second part

outpath<-"/Users/vborda/Documents/Public_data/INS_data/Chromopainter_runs/globetrotter/"
write.table(admixturesources_subset,file=paste0(outpath,"AncestryProportions_Sources_GLOBETROTTER_20251025.txt"),
            row.names=FALSE,quote=FALSE,col.names = TRUE)

SupplementaryData5<-pltable[,c("Cluster","date.1D","Date.CI","max.R1.","FQ1","FQ2")]
write.table(SupplementaryData5,file=paste0(outpath,"SupplementaryData5_GLOBETROTTER_20251025.txt"),
            row.names=FALSE,quote=FALSE,col.names = TRUE,sep = "\t")

###########################################################
## SORT OUT RESULTS TO PLOT
pltable <- final.res2plot
pltable <-pltable[popplot2,]
pltable$index <- as.numeric(row.names(pltable))       ### To Oder
pltable<-pltable[order(pltable$index), ]             ### To Oder

pltable <- pltable[!pltable$Result%in%c("U","NA"),]
pltable$Result[pltable$Result=="1D(2D)"] <- "1D"
pltable$Result[pltable$Result=="1MW(2D)"] <- "1MW"
#pltable<-na.omit(pltable) 

###################################
library(stringr)

rev_pops <- c("Chopccas","Matses")
rev_pops<-c()
pops <- as.character(pltable$Cluster)
noprop <- FALSE

tempresults <- addGTresults(pltable,rev_pops)
all_dates <- tempresults[[1]]
all_plot_mat <- tempresults[[2]]
all_src_mat <- tempresults[[3]]
dboots <- tempresults[[4]]
allsources <- tempresults[[5]]

names_pops<-row.names(all_dates)
names_order<- gsub("_a","",names_pops)

ii <-as.numeric(all_dates[,1])

makeGenFromDate <- function(x, year0=2010, gen_length=26.9){
  y <- round((((year0-x)/gen_length)-1),0)
  return(y)
}
##################################################################################
## PLOT

jpeg(filename = "/Users/vborda/Documents/Public_data/INS_data/Manuscript/Figures/Figure_2_Globetrotter_20251025.jpg",
     width = 20, height = 14, units = "cm", pointsize =12,  res = 600)

plot1<-layout(matrix(c(1,1,1,2,2,2,3,3,3),3),
              widths=c(4,2,2),heights=c(6,6,6))

topmar <- 6
n_pops <- 14

par(mar=c(4.5,12,topmar-1.5,0.5))
d_pops <- gsub("\\_a","",rownames(all_dates))
x_labs3<-c(3,6,9,12,15)
x_labs3char <-c("3","6","9","12","15")
x_max <- min(x_labs3)


poplabpos <- 1
ev1pos1 <- 2160
ev1pos2 <- 2140
ev2pos1 <- 2120
ev2pos2 <- 2100


## EMPTY PLOT FOR DATES
plot(0,0,xlim=range(x_labs3),  ###################     xlim=rev(range(x_labs3) con esto se volteaba el eje X
     ylim=c(n_pops,0.8),type="n",axes=F,xlab="",ylab="")   ###  ylim=c(n_pops,1.5) grosor del plot
axis(1,at=x_labs3,labels=x_labs3char,cex.axis=1,tick=F,padj=0.5)
for(j in x_labs3) abline(v=j,lty=2)
mtext("Date of Admixture (Generations)",1,line=3,cex = 0.8)

###################################################################
## PLOT DATES AND EVENT SOURCE ANCESTRY
## RUN THE CODE TWICE, THE FIRST TIME GETS THE ORDERING AND THE 
## SECOND TIME ACTUALLY PLOTS THE DATA
for(run in 1:2)
{
   if(run == 1)
   {
      poporder2 <- as.character(popplot2[1:14])
      popordertab <- c()
    }
    if(run == 2)
    {
        #poporder2 <- popordertab[,1]

    }
    
    for(i in 1:nrow(all_dates))
   {
        pop <- rownames(all_dates)[i]
        #pop <- all_dates$Cluster[i]
        poppos <- (1:length(poporder2))[poporder2%in%gsub("\\_a","",pop)]
        res <- all_dates[pop,5]
        #res <- all_dates[all_dates$Cluster==pop,4]
        d <- as.numeric(all_dates[pop,1])
        dh <- as.numeric(all_dates[pop,3])
        dl <- as.numeric(all_dates[pop,2])
        prop <- as.numeric(all_dates[pop,4])
  
        
        if(run == 2) points(d,poppos,pch=20, cex=3)################################################## plota la data
        if(run == 2) lines(x=c(dl,dh),y=c(poppos,poppos)) ########################################### Extremos de los intervalos de confianza
        if(res %in% c("1D","1MW"))
        {
          pop2 <- gsub("\\_a","",pop)
          anc1 <- pltable[pltable$Cluster==pop2,"best.source1"] ##### llama la primera fuente del primer evento
          anc2 <- pltable[pltable$Cluster==pop2,"best.source2"]
          if(run == 1 & pop2 %in% rev_pops)
          {
            addancs <- c(pop,anc2,anc1,dh)
          } else if(run == 1)
          {
            addancs <- c(pop,anc1,anc2,dh)
          }
          if(run == 1) popordertab <- rbind(popordertab,addancs)
          
          if(run == 2 & pop %in% rev_pops)
          {
            anc1 <- pltable[pltable$Cluster==pop2,"best.source2"]
            anc2 <- pltable[pltable$Cluster==pop2,"best.source1"]
          }
          anc3 <- pltable[pltable$Cluster==pop2,"best.source1.ev2"]
          anc4 <- pltable[pltable$Cluster==pop2,"best.source2.ev2"]
          pcol1 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Population==anc1])]
          pcol2 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Population==anc2])]
          pcol3 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Population==anc3])]
          pcol4 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Population==anc4])] ########### original en los 4 popkey$Ethnic_Group==anc4

        }
        if(res %in% c("2D"))
        {
          pop2 <- gsub("\\_a","",pop)
          anc1 <- pltable[pltable$Cluster==pop2,"best.source1.date1"]
          anc2 <- pltable[pltable$Cluster==pop2,"best.source2.date1"]
          if(pop2 %in% rev_pops)
          {
            addancs <- c(pop,anc2,anc1,dh)
          } else if(run == 1)
          {
            addancs <- c(pop,anc1,anc2,dh)
          }
          if(run == 1) popordertab <- rbind(popordertab,addancs)
          if(run == 2 & pop %in% rev_pops)
          {
            anc1 <- pltable[pltable$Cluster==pop2,"best.source2"]
            anc2 <- pltable[pltable$Cluster==pop2,"best.source1"]
          }
          anc3 <- pltable[pltable$Cluster==pop2,"best.source1.date2"]
          anc4 <- pltable[pltable$Cluster==pop2,"best.source2.date2"]
          pcol1 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Population==anc1])]
          pcol2 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Population==anc2])]
          pcol3 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Population==anc3])]
          pcol4 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Population==anc4])]
        }
    }
}
if(run == 1)
{
  popordertab <- popordertab[grep("_a",popordertab[,1],invert = T),]
  popordertab <- cbind(popordertab,
                       as.character(sapply(popordertab[,2],function(x){
                         popkey$RegionM[popkey$Population==x]})))
  popordertab <- cbind(popordertab,
                       as.character(sapply(popordertab[,3],function(x){
                         popkey$RegionM[popkey$Population==x]})))
  popordertab[,5] <- factor(popordertab[,5],levels=ancreg_list)
  popordertab[,6] <- factor(popordertab[,6],levels=ancreg_list)
  popordertab <- popordertab[order(popordertab[,6],popordertab[,5],
                                   round((1950-as.numeric(popordertab[,4]))/29)),]
}


###################################################################
## PLOT Y-AXIS NAMES AND COLOURS
y_ax_cols <- c()
for(i in poporder2)  #for(i in poporder) 
{
    y_ax_cols <- c(y_ax_cols,i)
}
for(i in 1:n_pops)
{
  
    neworder=c("Iquitos","Tumbes","Lambayeque","Trujillo","Huaraz","PEL","Lima","El Carmen","Arequipa","Ayacucho","Moquegua","Tacna","Cusco","Juliaca")
    axis(2,pos=poplabpos,at=i,labels=neworder[i],  ### MODIFICADO by NEGRI  labels=poporder[i]
    #axis(2,pos=poplabpos,at=i,labels=y_ax_cols[i],  ### MODIFICADO by NEGRI  labels=poporder[i]
         col.axis="black",las=2,tck=0,lwd=0,line=0,cex.axis=1.2) #line=-0.5 col.axis##############################################################################3
}

####################### To insert the type of result (one date, multiway-admixture, two dates, uncertain)##################3
for(i in 1:n_pops)
{resultG=y_ax_cols[i]
  resultG2 <-subset(final.res2plot$Result,final.res2plot$Cluster==resultG)
  
  axis(2,pos=poplabpos+700,at=i,labels=resultG2,  ### MODIFICADO by NEGRI  labels=poporder[i]
       col.axis="black",las=2,tck=0,lwd=0,line=-0.5,cex.axis=1.3) #line=-0.5 col.axis##############################################################################3
}

####################### Lines
for(i in seq(0.5,(n_pops+1),1)) abline(h=i,lty=3,lwd=1)  ### horizontal line in the date plot

###################################################################
## PLOT ADMIXTURE SOURCES
## FIRST EVENTS
color<-popkey$Colour[15:95]
plot_mat <- matrix(0,nrow=nrow(all_plot_mat),ncol=length(poporder2))
colnames(plot_mat) <- as.character(poporder2)
for(i in 1:ncol(plot_mat))
{
    pop2 <- as.character(poporder2[i])
    if(pop2 %in% colnames(all_plot_mat)) plot_mat[,pop2] <- unlist(all_plot_mat[,pop2])
}
#par(mar=c(9,0,topmar+0.8,0.8))
par(mar=c(9.2,0,topmar-0.2,0.8))
plot(0,0,xlim=c(0,1),
     #ylim=c(bp[length(bp)],bp[1]),type="n",axes=F,xlab="",ylab="")
     ylim=c(14,1),type="n",axes=F,xlab="",ylab="")
barplot(plot_mat[,1:n_pops],col=c(color,"white",color),
        yaxt="n",xlab="",beside=F,
        main="",cex.main=0.75,border=NA,horiz=T,
        xaxt="n",add=T)
text(y=-1,x=0.52,labels="Inferred haplotypic composition\nof the two admixing sources",xpd=T,srt=0,adj=0.5,cex=1.2) ## adj center align 0.5

###################################################################
## LEGEND
#par(mar=c(1,0,4,0))
rows<-c(grep(c("null.1"),rownames(admixturesources2)),grep(c("null.2"),rownames(admixturesources2)))
admixturesources5<-admixturesources2[rows,]
sources = admixturesources5[,colSums(admixturesources5) > 0.05]
leginfo2plot<-leginfo2order[leginfo2order$Population %in% colnames(sources),]

par(mar=c(5.5,2,3.5,0)) ##(bottom, left, top, right))
plot(0,0,axes=F,xlab="",ylab="",type="n")
legend_text <- leginfo2plot$Population
l <- legend("top",legend=leginfo2plot$Population, pch=c(rep(15,nrow(leginfo2plot))),
            col=leginfo2plot$Colour,bty="n",
            ncol=2,xpd=T,pt.cex=2,x.intersp=1.5,y.intersp=1.3,
            pt.lwd=3,cex=1, title="Donor Populations")

###################################################################

dev.off()



