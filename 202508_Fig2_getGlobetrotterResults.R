#########################################################################
### SCRIPT TO PULL OUT GLOBETROTTER RESULTS FROM MULTIPLE OUTPUTS     ###
#########################################################################
### Original Author: George Busby (https://github.com/georgebusby)
### Modified by :  Victor Borda - IGS - (January 2025)

### This script call other 4 R scripts : getLines.R, makeDate.R , getGlobetrotter.R and getGlobetrotter.R

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
library(roxygen2)
library(xtable)
library(devtools)

main_dir <-"/Users/vborda/Documents/Public_data/INS_data/Chromopainter_runs/globetrotter/globetrotter_Nov2024/Results/"
boot_dir <-"/Users/vborda/Documents/Public_data/INS_data/Chromopainter_runs/globetrotter/globetrotter_2025/"
GT_dir <-"/Users/vborda/Documents/general_pipelines/plotting_scripts/globetrotter_scripts/"

source(paste0(GT_dir,"getLines.R"))
source(paste0(GT_dir,"makeDate_Wang2023.R"))
source(paste0(GT_dir,"getGlobetrotter.R"))
#source(paste0(GT_dir,"getFastGlobetrotterDates.R"))
source(paste0(GT_dir,"getGlobetrotterDates_onefile.R"))
source(paste0(GT_dir,"convert.factors.to.strings.in.dataframe.R"))
###########################################################
## DEFINE VARIABLES AND FILES
res_dir <- main_dir #paste0(main_dir,"globetrotter_results/")           ## Results directory
popkey_file <- "/Users/vborda/Documents/Public_data/INS_data/Chromopainter_runs/globetrotter/PGP_PopulationKey_500K_August2025.txt"
filename <- ".main"                    			## Add the word of the file name pop_filename example: Awajun_500K_allD.null.txt , select "_500_allD"
popkey <- read.table(popkey_file,header=T)
popkey$Ethnic_Group <- toupper(popkey$Population)                 ## Add a column
popplot <- as.character(popkey$Population)      
dates_pops <-popkey$BirthDate

## creating a character variable of population names

##########################################################
##########################################################
## DEFINE VARIABLES TO STORE MAIN RESULTS
resmat <- resmatnull <- resmatalt <- matrix(NA,nrow=0,ncol=22)         	## Creating empty matrix for Results
dateboots <- matrix(NA,nrow=0,ncol=4)                                  	## Creating empty matrix for ONE DATE Bootstrap results
date2boots <- matrix(NA,nrow=0,ncol=4)                                 	## Creating empty matrix for TWO DATES Bootstrap results
admixturesources2 <- matrix(nrow=0,ncol=length(popplot))				        ## Creating Matrix to describe the specific contribution of a donor to a target population
colnames(admixturesources2) <- popplot
admixturesources3 <- admixturesources4 <- matrix(nrow=0,ncol=length(popplot)+20)	## Creating Matrix to describe parameters values and the specific contribution of a donor to a target population

##########################################################
## GET RESULTS
for(pop in popplot[1:14])		## Set the number (N) of populations to analyze. It is recommended that this number correspond to the first N populations of the list 
{
    print(paste("getting GLOBETROTTER results for:",pop))
    pop1 <- pop

    donors <- gsub("-",".",popplot)
    donors
    a <- "null"
    infile_null <- paste0(res_dir,pop,filename,".txt")
    infile_null
    infile_dates_null <- paste0(boot_dir,pop,".onedate.boot") ## just prefix
    infile_2dates_null <- paste0(boot_dir,pop,".2dates.boot")
    if(file.exists(infile_null))
    {
      ## GET NULL DATA
      ll <- getGlobetrotter(infile_null,donors)
      ## GET ONE DATE BOOTSTRAPS
      ld <- getGlobetrotterDates(infile_dates_null)
      ## GET TWO DATE BOOTSTRAPS
      ld2 <- getGlobetrotterDates(infile_2dates_null)
      
      if(nrow(ld)>0) dateboots <- rbind(dateboots,cbind(pop,"null",ld))
      if(nrow(ld2)>0) date2boots <- rbind(date2boots,cbind(pop,"null",ld2))
      ## ADD TO THE MAIN RESULTS TABLE
      resmat <- rbind(resmat,c(pop,ll[[1]],"null"))
      
    }
    
    ## PULL OUT THE SOURCES AND GENERATE COPYING VECTORS USING PREDMAT
    if(!is.null(ll[[2]]))
    {
      finn <- ll[[2]]
      ## REORDER finn TO MATCH PREDMAT; THERE ARE 8 SOURCES THAT TWO FOR EACH TYPE OF EVENT
      cvs <- finn[,15:ncol(finn)]   ### Considering that the list have 13 target populations
      rownames(cvs) <- paste(paste(pop,a,sep="."),1:8,sep=".")
      ## ADD TO THE MAIN MATRIX
      while(ncol(cvs)!=length(popplot)) cvs <- cbind(cvs,0)
      colns <- gsub("-",".",popplot)
      colnames(cvs)[colnames(cvs) == ""] <- colns[!colns%in%colnames(cvs)]
      admixturesources2 <- rbind(admixturesources2,cvs[,colns])
      ll_tmp <- c()
      for(k in 1:nrow(cvs))
      {
        ll_tmp <- rbind(ll_tmp,ll[[1]])  
      }
      admixturesources3 <- rbind(admixturesources3,cbind(ll_tmp,cvs[,colns]))
    }
 }

resmat <- data.frame(resmat,stringsAsFactors = FALSE)
resmat$result <- toupper(resmat$result)


#############################################################
## NOW WE HAVE OUR RESULTS TABLE AND WE CAN WORK WITH IT TO 
## GENERATE OUR FINAL RESULTS
#############################################################
## GEN ADMIXTURE P-VALUE BASED ON NULL BOOTSTRAPS
resmat$pval <- 1 ## modified original value 1
resmat$nump <- 0
nullruns <- grep("null",paste(dateboots$pop,dateboots[,2]))
nullruns <- dateboots[nullruns,]
for(i in as.character(unique(nullruns$pop)))
{
    bb <- nullruns$date1.est.boot[nullruns$pop==i]
    rowindex <- paste(resmat$V1,resmat$V22,sep=".") == paste(i,"main",sep=".")
    resmat$pval[rowindex] <- sum(bb>400 | bb < 1) / length(bb)
    resmat$nump[rowindex] <- length(bb)
    rowindex <- paste(resmat$V1,resmat$V22,sep=".") == paste(i,"null",sep=".")
    resmat$pval[rowindex] <- sum(bb>400 | bb < 1) / length(bb)
}


#############################################################
## GENERATE A FINAL TABLE OF RESULTS
## IMPORTANTE PARECE QUE LOS QUE SON DIFERENTE DE MULTIWAY NO LOS RECONOCE

##restypes <- c("MULTIPLE-DATES","ONE-DATE","ONE-DATE-MULTIWAY","NO-ADMIXTURE" ,"UNKNOWN") #editado

restypes <- c("ONE-DATE","ONE-DATE-MULTIWAY","MULTIPLE-DATES","UNCERTAIN" ,"NO-ADMIXTURE" ,"UNKNOWN")
numericcols <- c("maxR2fit.1date","fit.quality.1event","fit.quality.2events",
                 "maxScore.2events","proportion.date1.source1","proportion.date2.source1")
finaltable <- matrix(NA,nrow=0,ncol=ncol(resmat))
for(i in restypes)
{
  ii <- resmat[resmat$result==i,]
  if(nrow(ii) > 0)
  {
    ii[,numericcols] <- signif(apply(ii[,numericcols],2,as.numeric),3)
    #ii[,numericcols] <- apply(signif(apply(ii[,numericcols],2,as.numeric),3),2,as.character)
    ii <- convert.factors.to.strings.in.dataframe(ii)
    ii <- ii[order(ii$bestmatch.event1.source1,ii$gen.1date,decreasing=T),]
  }
  finaltable <- rbind(finaltable,ii)
}

finaltable <- cbind(finaltable[,c(1,22)],finaltable[,c(2:21,23:(ncol(finaltable)))])
colnames(finaltable) <- c("Cluster","Analysis","Result","Date","alpha",
                          "max(R1)","FQ1","FQ2","best.source1","best.source2",
                          "alpha2","best.source1.ev2","best.source2.ev2",
                          "Date.2a","Date2b","maxScore2ev","alpha2.date1",
                          "best.source1.date1","best.source2.date1",
                          "alpha2.date2","best.source1.date2","best.source2.date2","pval","nump")

#############################################################
## LOOK AT DATE BOOTSTRAPS AND GENERATE CONFIDENCE INTERVALS
dateboots2 <- dateboots3 <- c()
for(i in unique(popplot[1:14]))
{

    for(analy in c("null")) #original for(analy in c("main","null"))
    {
        ds <- round(dateboots$date1.est.boot[dateboots$pop==i & dateboots[,2] == analy])
        ### I am changing this line to keep the generation time
        #ii <- sapply(quantile(ds,c(0.975,0.025)),makeDate)
        ii <- round(quantile(ds,c(0.025,0.975)))
        dateboots2 <- rbind(dateboots2,c(i,analy,ii))
        ii <- round(quantile(ds,c(0.025,0.975)))
        dateboots3 <- rbind(dateboots3,c(i,analy,ii))
    }
}


#############################################################

### PARA DOS DATAS
### LOOK AT TWO DATE BOOTSTRAPS AND GENERATE CONFIDENCE INTERVALS
date2boots2 <- date2boots3 <- c()
for(i in unique(date2boots$pop))
{
    for(analy in c("null")) #original for(analy in c("main","null"))
    {
        ds <- round(date2boots$date1.est.boot[date2boots$pop==i & date2boots[,2] == analy])
        ds1 <- round(date2boots$date1.est.boot[date2boots$pop==i& date2boots[,2] == analy])
        ds2 <- round(date2boots$date2.est.boot[date2boots$pop==i& date2boots[,2] == analy])
        #ii <- sapply(quantile(ds1,c(0.975,0.025)),makeDate) ## original
        #iii <- sapply(quantile(ds2,c(0.975,0.025)),makeDate) ## original
        
        ii <- round(quantile(ds1,c(0.025,0.975)))  #I added 
        iii <- round(quantile(ds2,c(0.025,0.975))) #I added 
          
        date2boots2 <- rbind(date2boots2,c(i,analy,ii,iii))
        ii <- round(quantile(ds1,c(0.025,0.975)))
        iii <- round(quantile(ds2,c(0.025,0.975)))
        date2boots3 <- rbind(date2boots3,c(i,analy,ii,iii))
    }
}

#finaltable <- convert.factors.to.strings.in.dataframe(finaltable)  ORIGINAL
finaltable[] <- lapply(finaltable, as.character) ################### MODIFIED BY NEGRI

finaltable$Date.CI <- finaltable$dateL <- finaltable$dateH <- NA

#############################################################
## FILL IN DATES AND CIS IN THE MAIN RESULTS TABLE

rowindex <- match(paste(dateboots2[,1],dateboots2[,2],sep="."),
                      paste(finaltable$Cluster,finaltable$Analysis,sep="."),nomatch=0)
finaltable$Date.CI[rowindex] <- paste0("(",dateboots2[,3],"-",dateboots2[,4],")")
finaltable$dateL[rowindex] <- dateboots3[,3]
finaltable$dateH[rowindex] <- dateboots3[,4]

finaltable$Date2a.CI <- finaltable$date2aL <- finaltable$date2aH <- finaltable$Date2b.CI <- finaltable$date2bL <- finaltable$date2bH <- NA
rowindex <- match(paste(date2boots2[,1],date2boots2[,2],sep="."),
                  paste(finaltable$Cluster,finaltable$Analysis,sep="."),nomatch=0)

finaltable$Date2a.CI[rowindex] <- paste0("(",date2boots2[,3],"-",date2boots2[,4],")")
finaltable$date2aL[rowindex] <- date2boots3[,3]
finaltable$date2aH[rowindex] <- date2boots3[,4]
finaltable$Date2b.CI[rowindex] <- paste0("(",date2boots2[,5],"-",date2boots2[,6],")")
finaltable$date2bL[rowindex] <- date2boots3[,5]
finaltable$date2bH[rowindex] <- date2boots3[,6]

#############################################################
## ALSO ROUND DATE COLUMN
finaltable$Date <- round(as.numeric(finaltable$Date))
finaltable$Date.2a <- round(as.numeric(finaltable$Date.2a))
finaltable$Date2b <- round(as.numeric(finaltable$Date2b))
#############################################################
### MAKE SOME PRETTY RESULTS TABLES AND STORE FOR PLOTTING
allcols <- c("Cluster","Analysis","pval","maxScore2ev","max.R1.","FQ1","FQ2","Result")
onedatecols <- c(allcols,"Date","dateL","dateH","alpha","best.source1","best.source2")
onemultcols <- c(onedatecols,"alpha2","best.source1.ev2","best.source2.ev2") 
twodatecols <- c(allcols,"Date.2a","date2aH","date2aL","best.source1.date1","best.source2.date1","alpha2.date1","Date2b","date2bH","date2bL","best.source1.date2","best.source2.date2","alpha2.date2")

#############################################################
#############################################################
#############################################################
## NOW ALSO DEFINE NO-ADMIXTURE WHEN THE REDUCTION IN THE NULL
## INFERENCE R^2 IS GREATER THAN 1/3
###res.tabA <- convert.factors.to.strings.in.dataframe(finaltable) #[,twodatecols] ## ORIGINAL
res.tabA <- finaltable   ############################ NEGRI OPTION
res.tabA[] <- lapply(finaltable, as.character) ##### MODIFIED by NEGRI

finaltable[] <- lapply(finaltable, as.character)

################################################ NOT NECESSARY for NATIVES
## DEFINE NO ADMIXTURE FROM P-VALUES

#res.tabA <- convert.factors.to.strings.in.dataframe(finaltable) #[,twodatecols]
## DEFINE NO ADMIXTURE FROM P-VALUES
#res.tabA$Result[res.tabA$pval>0.05] <- "NO-ADMIXTURE"

##res.tabA$Result[res.tabA$pval>0.05] <- "NO-ADMIXTURE" ## FOR DIASPORA RESULT, I skip this line to show the non-significant results


res.tabA$Result[res.tabA$`max(R1)`>0.99] <- "ONE-DATE" ## Considering Hellenthal recommendations. If globetrotter classifies as 2D, if R1 is greater than 0.99 this should be considered one date

#############################################################
#############################################################
## ORDER RESULTS BY EVENT TYPE AND DATE
res.tabAll <- res.tabA
rt1 <- res.tabA[!res.tabA$Result%in%c("UNCERTAIN","NO-ADMIXTURE"),]
rt2 <- res.tabA[res.tabA$Result%in%c("UNCERTAIN","NO-ADMIXTURE"),]
res.ord <- rt1$best.source1.date1
rt1 <- rt1[order(res.ord,as.numeric(rt1$Date.2a)),]
res.ord <- c(res.ord,rep("XXX",nrow(rt2)))
res.tabA <- rbind(rt1,rt2)
res.res <- res.tabA$Result
res.res <- gsub("ONE-DATE-MULTIWAY","1MW",res.res)
res.res <- gsub("ONE-DATE","1D",res.res)
res.res <- gsub("MULTIPLE-DATES","2D",res.res)
res.res <- gsub("UNCERTAIN","U",res.res)
res.res <- gsub("NO-ADMIXTURE","NA",res.res)
res.res <- gsub("UNKNOWN","NA",res.res)
res.tabA$Result <- res.res
res.tabA[,"maxScore2ev"] <- round(as.numeric(as.character(res.tabA[,"maxScore2ev"])),2)
res.tabA[,"pval"] <- round(as.numeric(as.character(res.tabA[,"pval"])),2)


#############################################################
## SWITCH RESULT IF MULTIPLE DATE ARE UNREASONABLE ie CI IS LESS THAN 2
res.tabA$FinalResult <- res.tabA$Result
res.tabA$FinalAnaly <- "null"
min_gens <- 0
pops <- popplot
test1mw <- res.tabA$Cluster%in%pops & res.tabA$Result == "2D" & res.tabA$FQ1 <= 0.975
res.tabA$FinalResult[test1mw] <- "1MW"
test1mw <- res.tabA$Cluster%in%pops & res.tabA$Result == "2D" & res.tabA$FQ1 > 0.975
res.tabA$FinalResult[test1mw] <- "1D"
res.tabA$FinalAnaly[res.tabA$Cluster%in%pops] <- "null"

## This will change the columns
## SWAP SECOND DATE COLUMNS FOR FIRST DATE
#temp <- res.tabA[res.tabA$Cluster%in%pops,]
#date1cols <- c("Date","alpha","best.source1","best.source2","dateH","dateL","Date.CI")
#date2cols <- c("Date2b","alpha2.date2","best.source1.date2","best.source2.date2","date2bH","date2bL","Date2b.CI")
#res.tabA[res.tabA$Cluster%in%pops,date1cols] <- temp[,date2cols]
#res.tabA[res.tabA$Cluster%in%pops,date2cols] <- temp[,date1cols]

###################################################

## THINK ABOUT COMPARING TWO DATE INF TO ONE DATE
## WHAT ARE THE TWO DATES?
## HOW DOES THE SECOND DATE RELATE TO FIRST
## WHAT ARE THE SOURCES? ARE THEY THE SAME FOR BOTH EVENTS?
## IGNORE 
test <- res.tabA$Result=="2D"
tmp.res <- res.tabA[test==T,]
tmp.res <- tmp.res[grep("null",tmp.res$Analysis,invert=T),]
test <- (as.numeric(tmp.res$date2aL)>min_gens)
test[is.na(test)] <- F
tmp.res <- tmp.res[test==F,]
if(nrow(tmp.res)>0)
{
    tmp.res$FinalResult <- "1D"
    tmp.res$FinalResult[tmp.res$FQ1 <= 0.975] <- "1MW"
    res.tabA[rownames(res.tabA)%in%rownames(tmp.res),] <- tmp.res
}

## then generate a final table will all results
## plus the final results table 
colnames(res.tabA)[colnames(res.tabA)=="Date"] <- "date.1D"
colnames(res.tabA)[colnames(res.tabA)=="dateL"] <- "date.1D.L"
colnames(res.tabA)[colnames(res.tabA)=="dateH"] <- "date.1D.H"
colnames(res.tabA)[colnames(res.tabA)=="Date.2a"] <- "date.2D.1"
colnames(res.tabA)[colnames(res.tabA)=="date2aH"] <- "date.2D.1.H"
colnames(res.tabA)[colnames(res.tabA)=="date2aL"] <- "date.2D.1.L"
colnames(res.tabA)[colnames(res.tabA)=="Date2b"] <- "date.2D.2"
colnames(res.tabA)[colnames(res.tabA)=="date2bH"]  <- "date.2D.2.H"
colnames(res.tabA)[colnames(res.tabA)=="date2bL"] <- "date.2D.2.L"


#res.tabA <- res.tabA[order(res.tabA$FinalResult,res.tabA$best.source1),]
#AQUI
###
outputname<-"globetrotter_results_Admixed_PGP_20250812"

write.csv(res.tabA,file= paste0(main_dir,outputname,".csv"),quote=F,row.names=F)


#############################################################
## NOW WORK ON GETTING MAKING NICE TEX TABLES
res.tabB <- res.tabA
res.tabB$Result <- res.tabB$FinalResult
rt1 <- res.tabB[!res.tabB$FinalResult%in%c("U","NA"),]
rt2 <- res.tabB[res.tabB$FinalResult%in%c("U","NA"),]
res1.ord <- rt1[,"best.source1"]
res1.ord[rt1$Result=="2D"] <- rt1$best.source1.date1[rt1$Result=="2D"]
res.ord.levs <- popplot
res1.ord <- factor(res1.ord,levels=res.ord.levs)
res2.ord <- rt1[,"best.source2"]
res2.ord[rt1$Result=="2D"] <- rt1$best.source2.date1[rt1$Result=="2D"]
res2.ord <- factor(res2.ord,levels=res.ord.levs)

rt1 <- rt1[order(res2.ord,res1.ord,as.numeric(rt1[,"date.1D"])),]
res.ord <- c(res.ord,rep("XXX",nrow(rt2)))
res.tabB <- rbind(rt1,rt2)
## PULL OUT THE FINAL RESULTS
#rows2keep <- grep("main.null",res.tabB[,2],invert=F) ## OLD VERSION
res2keep <- res.tabB[,c("Cluster","FinalAnaly")]
res2keep <- res2keep[duplicated(res2keep),]
rows2keep <- c()
for(i in 1:nrow(res2keep))
{
    keeper <- (1:nrow(res.tabB))[res.tabB$Cluster==res2keep[i,1] & res.tabB$Analysis==res2keep[i,2]]
    rows2keep <- c(rows2keep,keeper)
}

#res.tabB <- res.tabB[rows2keep,]

#############################################################
## THE FINAL MAIN RESULTS TABLE
write.csv(res.tabB,file=paste0(main_dir,outputname,"_Alternative.csv"),quote=F)
#############################################################
res.tabcols <- c("Cluster","Analysis","pval","Result","FinalResult","max.R1.","FQ1","FQ2","maxScore2ev" ,
                 "date.1D","Date.CI", "alpha","best.source1","best.source2",
                 "alpha2","best.source1.ev2","best.source2.ev2",
                 "date.2D.1","Date2a.CI","alpha2.date1","best.source1.date1","best.source2.date1",
                 "date.2D.2","Date2b.CI","alpha2.date2","best.source1.date2","best.source2.date2")
colnames(res.tabB)[colnames(res.tabB)=="max(R1)"] <- "max.R1."
res.tabB <- res.tabB[,res.tabcols]

test <- res.tabB$Result!=res.tabB$FinalResult
res.tabB$Result[test] <- paste(res.tabB$FinalResult[test],"(",res.tabB$Result[test],")",sep="")
res.tabB$Result[test] <- res.tabB$FinalResult[test]
res.tabB <- res.tabB[,colnames(res.tabB)!="FinalResult"]

for(i in c("date.1D","date.2D.1","date.2D.2"))
{
    res.tabB[,i] <- sapply(res.tabB[,i],function(x) {makeDate(round(as.numeric(x)),add_BCE=F)})
    res.tabB[,i][which(res.tabB[,i]<0)] <- paste0(-res.tabB[,i][which(res.tabB[,i]<0)],"B")
    res.tabB[,i] <- as.character(res.tabB[,i])
}

#############################################################
## SOME FURTHER COSMETIC CHANGES:: change p = 0 to p<0.01
res.tabB$pval[res.tabB$pval==0] <- "<0.01"
## only include best guess info
onedatecols <- c("date.1D","Date.CI","alpha","best.source1","best.source2" )
onemwcols <- c( "alpha2","best.source1.ev2","best.source2.ev2")
twodatecols <- c("date.2D.1","Date2a.CI","alpha2.date1","best.source1.date1","best.source2.date1",
                 "date.2D.2","Date2b.CI","alpha2.date2","best.source1.date2","best.source2.date2")

res.tabB[res.tabB$Result%in%c("1D","1D(2D)"),c(onemwcols,twodatecols)] <- NA
res.tabB[res.tabB$Result%in%c("1MW","1MW(2D)"),twodatecols] <- NA
res.tabB[res.tabB$Result=="2D",c(onedatecols,onemwcols)] <- NA

final.res2plot <- res.tabB

write.table(res.tabA,file=paste0(main_dir,outputname,"_toplot.txt"))
write.table(admixturesources2,paste0(main_dir,outputname,"_AdmixtureSources2.txt"))
write.table(admixturesources3,paste0(main_dir,outputname,"_AdmixtureSources3.txt"))
write.table(admixturesources4,paste0(main_dir,outputname,"_AdmixtureSources4.txt"))
write.table(dateboots,paste0(main_dir,outputname,"_OneDateBootstraps.txt"))
write.table(date2boots,paste0(main_dir,outputname,"_TwoDateBootstraps.txt"))


