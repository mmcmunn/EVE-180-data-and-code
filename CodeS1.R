#this script is an analysis of UC Davis EVE 180 class data from 2017
#it contains all final analyses for the project
#Marshall McMunn is the primary author of the script: mmcmunn@gmail.com
#clear all objects
#rm(list=ls())

#LOAD PACKAGES
#####
library("tidyverse")
library("lme4")
library("vegan")
library("reshape2")
library("mvabund")
library("RColorBrewer")
library("gridExtra")
library("betareg")
library("lmtest")
library("MASS")
#####
#standard error function
se <- function(x) sqrt( var(x[ !is.na(x) ] ) / length(x[ !is.na(x) ] ))
#LOADING DATA
getwd()
#####
#set working directory to same location as this script file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load data
  #a night survey of spider abundance on cages, twice in riparian block, once in other blocks
    nightSpiD <- read.csv("DataS1.csv",stringsAsFactors = FALSE)

    
  #all plant data - heribivory, size, dates, species, twice weekly
    plantD <- read.csv("DataS2.csv",stringsAsFactors = FALSE)
    
  #a single night predation experiment performed across all replicates of fly disappearance rate during night and day
    predExpD <- read.csv("DataS3.csv",stringsAsFactors = FALSE)
  
  # s = 24 hour sticky traps, performed weekly at each replicate
  # p = 24 hour pitfall traps, performed weekly, pitfall collections by replicate
    insectD <- read.csv("DataS4.csv",stringsAsFactors = FALSE)

  #a list of all treatments, confirmed typo free for reference and consistency across datasets
    trtMaster <- read.csv("DataS5.csv", stringsAsFactors = FALSE)
    
  #reference for predator non-predtor classification
    predNonPred <- read.csv("DataS6.csv", stringsAsFactors = FALSE)
#####


#DATA FORMATTING

#lump LP and L from all datasets
    #LP represents lit predator removal. This treatment failed to reduce
    #predators overall, band did not alter overall arthropod abundance   
    plantD[plantD$trt=="LP" , "trt" ] <- "L"
    nightSpiD[nightSpiD$Treatment=="LP" , "Treatment" ] <- "L"
    predExpD[predExpD$Treatment=="LP" , "Treatment"] <- "L"
    insectD[insectD$treatment=="LP" , "treatment"] <- "L"
    trtMaster[trtMaster$treatment=="LP" , ] <- "L"
    
    
    
    
    
#####
#data cleaning for plant data
    plantD$propDam <- plantD$numDam/plantD$numLeaf
    plantD$species <- as.factor(plantD$species)
    plantD$trt <- as.factor(plantD$trt)
#add column to specify if measurement is initial or not
    plantD$initial <- ifelse(plantD$date=="2017-04-11" |plantD$date== "2017-04-25" |plantD$date== "2017-05-09",
           "init", "noninit")

    #add an ID column to keep track of each plant
    plantD$plantID<-paste(plantD$species , plantD$rep ,plantD$block,plantD$cohort , sep = ".")
    
#####
#ANALYSIS
#####
#percent area damage definitely not normal!
#check data
      table(plantD$date, plantD$cohort, plantD$alive)
#removing third week of data for cohort 1, this makes all measurements 1 week after planting
#for herbivory analysis, initial measuremesnts are removed here as well
      d <- plantD[plantD$initial=="noninit" & plantD$date!= "2017-04-24"&plantD$date!= "2017-05-02" , ]
      d<-d[!is.na(d$propDam),]
      
      tapply(d$propDam , list(d$trt , d$species, d$cohort), mean, na.rm = TRUE)
      tapply(d$propDam , list(d$trt), mean, na.rm = TRUE)
      tapply(d$areaDam , list(d$trt), mean, na.rm = TRUE)

#label treatments unlint and lit
      d$trt<-factor(d$trt , levels = c("U","L"))
      
      median(plantD[plantD$initial=="init" , "height"])
      sd(plantD[plantD$initial=="init" , "height"])
#prepare herbivory data for beta regresson
      dBeta <- d
  # calculate beta correction for sample size 447
      (1*(447-1) + 0.5)/447
      (0*(447-1) + 0.5)/447
      
     
      dim(dBeta)
      dBeta[dBeta$propDam==1 , "propDam"] <- 0.999
      dBeta[dBeta$propDam==0 , "propDam"] <- 0.001
      ?betareg
#proportion leaves damaged
      M1<-betareg(propDam ~ trt * species + block + as.factor(cohort), data = dBeta)
      summary(M1)
      
      M1.1<-betareg(propDam ~  trt + species + block + as.factor(cohort), data = dBeta)
      lrtest(M1, M1.1)
      summary(M1.1)
      
      M1.2<-betareg(propDam ~  species + block + as.factor(cohort), data = dBeta)
      lrtest(M1.1, M1.2)
#separate analyses by species   
      #tomato
            unique(d$species)
            M1<-betareg(propDam ~ trt + block + as.factor(cohort), data = dBeta[dBeta$species=="t",])
            summary(M1)
            
            M1.2<-betareg(propDam ~ block + as.factor(cohort), data = dBeta[dBeta$species=="t",])
            lrtest(M1, M1.2)
            
            
      #brassica
            unique(d$species)
            M1<-betareg(propDam ~ trt + block + as.factor(cohort), data = dBeta[dBeta$species=="b",])
            summary(M1)
            
            M1.2<-betareg(propDam ~ block + as.factor(cohort), data = dBeta[dBeta$species=="b",])
            lrtest(M1, M1.2)     
            
          
      #pea
            unique(d$species)
            M1<-betareg(propDam ~ trt + block + as.factor(cohort), data = dBeta[dBeta$species=="p",])
            summary(M1)
            
            M1.2<-betareg(propDam ~ block + as.factor(cohort), data = dBeta[dBeta$species=="p",])
            lrtest(M1, M1.2)    
 
  #proportion area damaged
  #prepare area for beta regression
      dBeta$areaDam <- d$areaDam/100
      dBeta <- dBeta[- which(is.na(dBeta$areaDam)) , ]
      dBeta[dBeta$areaDam==1 , "areaDam"]<- 0.999
      dBeta[dBeta$areaDam==0 , "areaDam"] <- 0.001
      
      M1<-betareg(areaDam ~ trt * species + block + as.factor(cohort), data = dBeta)
      summary(M1)
      
      M1.1<-betareg(areaDam ~  trt + species + block + as.factor(cohort), data = dBeta)
      lrtest(M1, M1.1)
      summary(M1.1)
      
      M1.2<-betareg(areaDam ~  species + block + as.factor(cohort), data = dBeta)
      lrtest(M1.1, M1.2)

      
      #no difference in plant height
      M2 <- with(d , lm(height ~ trt + block + species + cohort))
      M2.1 <- with(d , lm(height ~  block + species + cohort))
      summary(M2)
      anova(M2, M2.1)
      
      #no difference in plant leaf number
      M3 <- with(d , glm.nb(numLeaf ~ trt + block + species + cohort))
      M3.1 <- with(d , glm.nb(numLeaf ~ block + species + cohort))
      
      summary(M3)
      anova(M3, M3.1)
      
      #no difference in plant leaf area
      M4 <- with(d , lm(areaLeaf ~ trt + block + species + cohort))
      M4.1 <- with(d , lm(areaLeaf ~ block + species + cohort))
      
      summary(M4) 
      anova(M4, M4.1)
      
      #no difference in survival
      M5 <- with(d , glm(alive ~ trt + block + species + cohort) , family = "binomial")
      M5.1 <- with(d , glm(alive ~  block + species + cohort) , family = "binomial")
      summary(M5)

      
#for plotting, add rows for overall data as a "species" allowing gg faceting to include
      #sumamry and each species side-by-side
    plantD_plotting <- plantD
    plantD_plotting[ , "species"] <- "Overall"
    plantD_plotting <- rbind(plantD_plotting , plantD)
    plantD_plotting$species <- factor(plantD_plotting$species , levels = c("Overall", "b" , "p", "t"))
    levels(plantD_plotting$species) <- c("Overall" , "Brassica" , "Pea" , "Tomato")
      
# Figure 4 - proportion Damage - by species
      P1 <- ggplot(plantD_plotting , aes(x = trt, y = propDam)) + facet_wrap(~species, ncol = 4)
      P1<-P1 + geom_bar(stat = "summary", fun.y = "mean", 
                    fill = c("#E69F00", "#999999", "#E69F00", 
                             "#999999", "#E69F00", "#999999",
                             "#E69F00", "#999999")) +
        stat_summary(fun.data = mean_se,geom = "errorbar", width = 0.2) +
        labs(x = "treatment", y = "proportion leaves damaged") +theme_bw() + scale_x_discrete(labels = c("ALAN" , "control"))

      pdf(file = "figures/Fig 4 - herbivory.pdf")
      P1
      dev.off()

#####
#ARTHROPOD SUMMARY STATISTICS AND PLOTS
####
#total abundance (sticky + pitfall)
      #how many records?
      dim(insectD)
      
      insectD$sampleID <- paste(insectD$block , insectD$replicate, insectD$date , sep = ".")
      abundTable <- data.frame(counts = as.vector(tapply(insectD$sampleID,insectD$sampleID, length)))
      
      abundTable$treatment <- insectD[ match(rownames(tapply(insectD$sampleID,insectD$sampleID, length)), insectD$sampleID) , "treatment"]
      abundTable$block <- insectD[ match(rownames(tapply(insectD$sampleID,insectD$sampleID, length)), insectD$sampleID) , "block"]
      #create a cohort column in insect dataframe
      insectD$cohort <- ifelse(insectD$date=="4/20/2017" , 1 , 0 )
      insectD$cohort <- ifelse(insectD$date=="4/26/2017" | insectD$date=="4/27/2017" | insectD$date=="5/4/2017" , 2 , insectD$cohort )
      insectD$cohort <- ifelse(insectD$date=="5/11/2017" , 3 , insectD$cohort )
      
  
      abundTable$cohort <- insectD[ match(rownames(tapply(insectD$sampleID,insectD$sampleID, length)), insectD$sampleID) , "cohort"]
      abundTable$s_p <- insectD[ match(rownames(tapply(insectD$sampleID,insectD$sampleID, length)), insectD$sampleID) , "s_p"]
#there was twice as much sampling during cohort 2 so these numbers are divided by two and rounded up to nearest integer
      abundTable$count <- ifelse(abundTable$cohort==2 ,ceiling(abundTable$count / 2) , abundTable$count)
  
      #overall arhtropod abundance models
      abundM0 <- with(abundTable , glm.nb(counts ~ treatment * block + treatment * as.factor(cohort) + treatment*s_p))
      abundM0.1 <- with(abundTable , glm.nb(counts ~ treatment + block + treatment * as.factor(cohort) + treatment*s_p))
      abundM0.2 <- with(abundTable , glm.nb(counts ~ treatment * block + treatment + as.factor(cohort) + treatment*s_p))
      abundM0.3 <- with(abundTable , glm.nb(counts ~ treatment * block + treatment * as.factor(cohort) + treatment+s_p))
      abundM0.4 <- with(abundTable , glm.nb(counts ~ treatment + block + treatment + as.factor(cohort) + treatment+s_p))
      abundM1.3 <- with(abundTable , glm.nb(counts ~block + as.factor(cohort) +s_p))
      
      #treatment X block
      anova(abundM0, abundM0.1, test = "Chisq")
      100* with(abundTable , tapply(counts , list(treatment , block) , mean ))[1,] /
        with(abundTable , tapply(counts , list(treatment , block) , mean ))[2,]  
      
      #treatment X cohort
      anova(abundM0, abundM0.2, test = "Chisq")
      100* with(abundTable , tapply(counts , list(treatment , cohort) , mean ))[1,] /
        with(abundTable , tapply(counts , list(treatment , cohort) , mean ))[2,]  
      
      #treatment X trap type
      anova(abundM0, abundM0.3, test = "Chisq")
     100* with(abundTable , tapply(counts , list(treatment , s_p) , mean ))[1,] /
        with(abundTable , tapply(counts , list(treatment , s_p) , mean ))[2,]
    
     
 
     #full model
      summary(abundM0)
      anova(abundM0.4, abundM1.3, test = "Chisq")
      
      100* with(abundTable , tapply(counts , list(treatment) , mean ))[1] /
        with(abundTable , tapply(counts , list(treatment) , mean ))[2]
      
#plot mean abundance and se - Fig 3a
      P1 <- ggplot(abundTable , aes(x = treatment, y = counts))
      P1 <-P1 + geom_bar(stat = "summary", fun.y = "mean", fill = c("#E69F00", "#999999")) +
        stat_summary(fun.data = mean_se,geom = "errorbar", width = 0.2) +
        labs(x = "treatment", y = "mean arthropods collected")  + scale_x_discrete(labels = c("ALAN" , "control")) +
        theme_bw()

      head(insectD)
      
#diversity in arthropod samples
      divTable <-  with(insectD, table(sampleID , order))
      divTable

      divSumTable <- data.frame(shannon = diversity(divTable) , 
                                alpha = specnumber(divTable),
                                evenness =  diversity(divTable)/log(specnumber(divTable)))
                                
      divSumTable$treatment <- insectD[ match(rownames(divSumTable), insectD$sampleID) , "treatment"]
      divSumTable$block <- insectD[ match(rownames(divSumTable), insectD$sampleID) , "block"]
      divSumTable$block <- insectD[ match(rownames(divSumTable), insectD$sampleID) , "cohort"]
      
      #shannon diversity the same
      anova(lm(shannon ~ treatment + block + treatment , data = divSumTable))
      
      #alpha diversity lower in unlit
      anova(lm(alpha ~ treatment + block + treatment, data = divSumTable))

      #evenness lower in unlit
      anova(lm(evenness ~ treatment + block + treatment , data = divSumTable))
      
      
      
      #plot means and se Shannon - Fig 3c
      P3 <- ggplot(divSumTable , aes(x = treatment, y = shannon))
      P3 <- P3 + geom_bar(stat = "summary", fun.y = "mean", fill = c("#E69F00", "#999999")) +
        stat_summary(fun.data = mean_se,geom = "errorbar", width = 0.2) +
        labs(x = "treatment", y = "mean Shannon diversity")  + scale_x_discrete(labels = c("ALAN" , "control")) +
        theme_bw()

      #plot means se of alpha diversity - Fig 3b
      P2 <- ggplot(divSumTable , aes(x = treatment, y = alpha))
      P2 <- P2 + geom_bar(stat = "summary", fun.y = "mean", fill = c( "#E69F00", "#999999")) +
        stat_summary(fun.data = mean_se,geom = "errorbar", width = 0.2) +
        labs(x = "treatment", y = "mean alpha diversity")  + scale_x_discrete(labels = c("ALAN" , "control")) +
        theme_bw()

      abundTable <- data.frame(tapply(insectD$sampleID, list(insectD$sampleID, insectD$order), length))
      abundTable[is.na(abundTable)] <- 0
      
      
#Supp Figure 3
      insectLong <- melt(divTable)
      sampleCodes<-strsplit(as.character(insectLong$sampleID) , split = "\\." )
      insectLong$date<-sapply(sampleCodes , "[", 3 )         
      insectLong[insectLong$date=="4/26/2017" , "date"] <- "4/27/2017"
      insectLong$trt<-insectD[match(insectLong$sampleID , insectD$sampleID) , c( "treatment" )]
      insectLong$date<-as.POSIXct(insectLong$date, format = "%m/%d/%Y")
    
      
     insectLong[insectLong$trt=="L" , "trt"] <- "ALAN"
     insectLong[insectLong$trt=="U" , "trt"] <- "Control"
     sum(insectLong$value)
      
      #plot - split by order
      pdf("../figures/suppFig2_time2.pdf" , width = 14, height = 12)
      
      P6<-ggplot(data = insectLong,
                 aes(x = date, y = value, color = trt)) +
        stat_summary(fun.y = mean, geom = "line")+stat_summary(fun.y = mean, geom = "point")+
      facet_wrap(vars(order), scales = "free")+
        ylab("mean abundance") + theme_bw() + scale_color_manual(values = c( "#E69F00", "#999999")) +
        guides(color=guide_legend(title="treatment"))
 P6
 
      dev.off()
      
      
     sort( table(insectD$order))
      
      
#proportion predator
      predTable <- as.data.frame.matrix(with(insectD, table(sampleID, predacious)))
      predTable$propPred <- predTable$Y/(predTable$N+predTable$Y)
      
      predTable$treatment <- insectD[ match(rownames(predTable), insectD$sampleID) , "treatment"]
      predTable$block <- insectD[ match(rownames(predTable), insectD$sampleID) , "block"]
      predTable$cohort <- insectD[ match(rownames(predTable), insectD$sampleID) , "cohort"]
      predTable$replicate <- insectD[ match(rownames(predTable), insectD$sampleID) , "replicate"]
      
      predTable[predTable$propPred==0 , "propPred"] <- 0.001
      predTable[predTable$propPred==1 , "propPred"] <- 0.999
      
      M2.1 <-betareg(predTable$propPred ~ predTable$treatment + predTable$block + predTable$cohort)
      summary(M2.1)
      head(predTable)
      with(predTable , tapply(propPred ,treatment  , mean))[1]/with(predTable , tapply(propPred ,treatment  , mean))[2]
      
      
      
      #predator abundance - yes more in lit
      M3.1<- glm.nb(Y~treatment, data = predTable)
      M3.0<- glm.nb(Y~1, data = predTable)
      anova(M3.1, M3.0, test = "Chisq")
      
      #non-predator abundance - yes, more in lit
      M4.1<- glm.nb(N~treatment, data = predTable)
      M4.0<- glm.nb(N~1,  data = predTable)
      anova(M4.1, M4.0, test = "Chisq")
      
     
      
      #proportion predator
      P4 <- ggplot(predTable , aes(x = treatment, y = propPred))
      P4 <- P4 + geom_bar(stat = "summary", fun.y = "mean", fill = c("#E69F00", "#999999")) +
        stat_summary(fun.data = mean_se,geom = "errorbar", width = 0.2) +
        labs(x = "treatment", y = "proportion predator") + scale_x_discrete(labels = c("ALAN" , "control")) +
        theme_bw()

      
      #community composition
      abundTable
      set.seed(12435324)
      ord.all <- metaMDS(comm=as.matrix(abundTable),distance="bray",
                         try = 100, k=2, noshare=(engine="isoMDS"))
      
      xy.coord <- as.data.frame(scores(ord.all))
      
      xy.coord$treatment <- insectD[ match(rownames(xy.coord) , insectD$sampleID) ,"treatment"  ]
      xy.coord$block <- insectD[ match(rownames(xy.coord) , insectD$sampleID) ,"block"  ]
      
      P5 <- ggplot(xy.coord, aes(x = NMDS1 , y = NMDS2, color = factor(treatment, labels = c("ALAN","control")))) +
        geom_point(size = 0.45)+scale_color_manual(values = c("#E69F00", "#999999")) +theme_bw()+ labs(color = "treatment")
      
      
    
      #yes, both are significant, and their interaction
      adonis(abundTable$counts ~ abundTable$treatment * abundTable$block + as.factor(abundTable$cohort) * abundTable$treatment)
      
      
#figure 2 - overall arthropods
      pdf("../figures/Fig 2 - overall arthropods.pdf")
      grid.arrange(P1 , P2 , P3 , P4, P5)
      dev.off()

        
#spider web abundance - takes some data wrangling because two surveys - one complete, one partial
      
      nightSpiD$obsLevelID <- paste(nightSpiD$Block ,nightSpiD$Replicate, nightSpiD$Date, sep = "." )
      spiderD <- data.frame(spiderCount = with(nightSpiD , tapply( Number.of.individuals,obsLevelID , sum)))
      spiderD$block <- nightSpiD[match(rownames(spiderD),nightSpiD$obsLevelID) , "Block"]
      spiderD$treatment <- nightSpiD[match(rownames(spiderD),nightSpiD$obsLevelID) , "Treatment"]
      spiderD$date <- nightSpiD[match(rownames(spiderD),nightSpiD$obsLevelID) , "Date"]
      spiderD$replicate <- nightSpiD[match(rownames(spiderD),nightSpiD$obsLevelID) , "Replicate"]
      
      
      spiderD[is.na(spiderD$spiderCount) , "spiderCount"] <- 0
      
      
#need to add zeros for unrecorded replicates - data entered were counts and didnt write down zeros
      spiderD511<- spiderD[spiderD$date=="5/11" , ]
      
      spiderD<- spiderD[spiderD$date=="5/17" , ]
      
      fullset <- data.frame(block = c(rep("R", 30), rep("B", 30), rep("G", 30)))
      fullset$replicate <- rep(1:30, 3)
      notRecorded <- fullset[!paste(fullset$block , fullset$replicate, sep = ".") %in% paste(spiderD$block , spiderD$replicate, sep = ".") , ]
      notRecorded$spiderCount <- 0
      notRecorded$date <- "5/17"
      notRecorded$treatment <-  insectD[ match( paste(notRecorded$block, notRecorded$replicate, sep = "."),
                                                paste(insectD$block , insectD$replicate, sep = ".")
      ) , "treatment" ]
      
      fullset511 <- data.frame(block = c(rep("R", 30)))
      fullset511$replicate <- 1:30
      notRecorded511 <- fullset511[!paste(fullset511$block , fullset511$replicate, sep = ".") %in% paste(spiderD511$block , spiderD511$replicate, sep = ".") , ]
      notRecorded511<-notRecorded511[notRecorded511$block=="R" , ]
      notRecorded511$spiderCount <- 0
      notRecorded511$date <- "5/11"
      notRecorded511$treatment <-  insectD[ match( paste(notRecorded511$block, notRecorded511$replicate, sep = "."),
                                                paste(insectD$block , insectD$replicate, sep = ".")
      ) , "treatment" ]
      
      spiderD <- rbind(notRecorded, spiderD, spiderD511, notRecorded511)
      
      #poisson model, note, data cannot be fit by negative binomial
      M6 <- glm(spiderCount~as.factor(treatment)*as.factor(block) ,family = "poisson", data = spiderD[spiderD$date=="5/17",])
      M6.1 <- glm(spiderCount~as.factor(treatment)+as.factor(block) ,family = "poisson", data = spiderD[spiderD$date=="5/17",])
      M6.2 <- glm(spiderCount~as.factor(treatment) ,family = "poisson", data = spiderD[spiderD$date=="5/17",])
      M6.3 <- glm(spiderCount~as.factor(block) ,family = "poisson", data = spiderD[spiderD$date=="5/17",])
      table(spiderD$date , spiderD$treatment)
      
      summary(M6.1)
      anova(M6, M6.1, test = "Chisq")
      anova(M6.2, M6.1, test = "Chisq")
      anova(M6.3, M6.1, test = "Chisq")
      
      #binomial model - spiders Y/N
      spiderD$present <- ifelse(spiderD$spiderCount==0 , 0 , 1)
      
      Mb1  <- glm(present~treatment+block,family = binomial, data = spiderD)
      Mb0  <- glm(present~block,family = binomial, data = spiderD)
      anova(Mb1, Mb0, test = "Chisq")
      with(spiderD[spiderD$date=="5/17",], tapply(present , treatment , mean))

      #plot means se of spider web data
      P1 <- ggplot(spiderD , aes(x = treatment, y = present))
      P1 + geom_bar(stat = "summary", fun.y = "mean", fill = c("#E69F00", "#999999")) +
        theme(axis.text.x = element_text(face = "bold", color = "black", size = 16)) +
        theme(axis.text.y = element_text(face = "bold", color = "black", size = 16)) +
        theme(axis.title = element_text(face = "bold", size = 18)) +
        labs(x = "treatment", y = "proportion with webs")

      
# response in abundance by Order - phototaxis
      abundTable <- data.frame(tapply(insectD$sampleID, list(insectD$sampleID, insectD$order), length))
      abundTable[is.na(abundTable)] <- 0
      
      treatments <- insectD[ match(rownames(abundTable) , insectD$sampleID) , "treatment"]
      cohorts <- insectD[ match(rownames(abundTable) , insectD$sampleID) , "cohort"]
      abundTable <- mvabund(abundTable)
      
      M5 <- manyglm(abundTable ~ as.factor(treatments), family = "negative.binomial")
      summary.manyglm(M5)
      anova.manyglm(M5)
      M5.sig <- anova(M5, p.uni = "adjusted")
      
      
    
      taxaRespond <- names(which(M5.sig$uni.p[2,]<=.05))
      summarizeOrder <- function(order, group = insectD$order){
        taxonSumm <- list()
        temp <- insectD[insectD$order == order, ]
        
        taxonSumm[["total abundance"]] <- nrow(temp)
        taxonSumm[["predacious"]] <- temp$predacious[1]
        taxonSumm[["ratio of pos. phototaxis"]] <- format( (nrow(temp[temp$treatment=="L", ])/2) /
          nrow(temp[temp$treatment=="U", ]), nsmall = 2)
        taxonSumm[["glm.nb p < "]] <-   M5.sig$uni.p[2,order]
        taxonSumm[["glm.nb (theta)"]] <- format(round(M5$theta[order] , 2) , digits = 2)
        taxonSumm[["glm.nb (2 * log likelihood)"]] <- format(round(M5$two.loglike[order] , 2) , digits = 2)
        taxonSumm
      }
      
      insectD$order <- gsub( " - " , "..." , insectD$order)
      
      orderSummary <- data.frame(sapply(unique(insectD$order) , function(x) unlist(summarizeOrder(x) )))
      orderSummary   <- orderSummary[ , order(-as.numeric(t(orderSummary[1,]) ) ) ]
      
      write.csv(print(orderSummary), file = "orderSummary.csv")
      
#
      
      orderMeans  <-  apply( abundTable , 2, function(x){ tapply( x, predTable$treatment  , mean )} )
      orderSEs <- apply( abundTable , 2, function(x){ tapply( x, predTable$treatment  , se )} )
      
      orderMeans<-melt(orderMeans)
      orderSEs<-melt(orderSEs)
      orderSumm <- orderMeans 
      
      orderSumm$SE <- orderSEs$value
      colnames(orderSumm) <- c("treatment", "order", "mean", "SE")

      orderSumm <- orderSumm[  orderSumm$order %in% taxaRespond , ]
      
      spiderSumm<-data.frame(cbind(
      with(spiderD , tapply(spiderCount ,treatment, mean)) ,
      with(spiderD , tapply(spiderCount ,treatment, se))
      ))
      spiderSumm<-cbind(spiderSumm , "Aranea - web census")
      spiderSumm <- cbind(spiderSumm , c("L" , "U"))
      spiderSumm<-spiderSumm[ , c(4, 3, 1, 2)]
      colnames(spiderSumm) <- colnames(orderSumm)
      
      
      
      orderSumm <- rbind(spiderSumm , orderSumm )
      pdf(file = "../figures/Fig 3 - sig order abundance.pdf",width = 10, height = 6)
      P1 <- ggplot(orderSumm , aes(x = treatment, y = mean, fill = treatment)) + facet_wrap(~order, ncol=4,scales="free")
      P1 + geom_bar(stat = "identity", fill = rep(c("#E69F00", "#999999"),11)) +scale_x_discrete(labels = c("ALAN" , "control"))+
        geom_errorbar(mapping = aes(x=treatment, ymin = mean-SE , ymax = mean + SE), width = 0.2)+
        labs(x = "treatment", y = "mean count per sample")
      dev.off()
      
      
      
      M6 <- manyglm(abundTable ~ as.factor(cohorts), family = "negative.binomial")
      summary.manyglm(M6)
      anova.manyglm(M6)
      M6.sig <- anova(M6, p.uni = "adjusted")
      
      taxaRespond <- names(which(M5.sig$uni.p[2,]<=.05))
      summarizeOrder <- function(order, group = insectD$order){
        taxonSumm <- list()
        temp <- insectD[insectD$order == order, ]
        
        taxonSumm[["total abundance"]] <- nrow(temp)
        taxonSumm[["predacious"]] <- temp$predacious[1]
        taxonSumm[["ratio of pos. phototaxis"]] <- format( (nrow(temp[temp$treatment=="L", ])/2) /
                                                             nrow(temp[temp$treatment=="U", ]), nsmall = 2)
        taxonSumm[["glm.nb p < "]] <-   M5.sig$uni.p[2,order]
        taxonSumm[["glm.nb (theta)"]] <- format(round(M5$theta[order] , 2) , digits = 2)
        taxonSumm[["glm.nb (2 * log likelihood)"]] <- format(round(M5$two.loglike[order] , 2) , digits = 2)
        taxonSumm
      }
      
      
      #supp figure 2 - method X order X predator
      insectDplot <- insectD
      insectDplot$s_p <- ifelse(insectDplot$s_p=="S" ,"Sticky Trap" , "Pitfall Trap" )
      png(file = "abundance X trap X order X pred.png", width = 10 , height = 5, units = "in", res = 300)
      P1 <- ggplot(insectDplot , aes(x = order, y = length(order), fill = predacious)) + facet_wrap(~s_p,scales="free") +
      scale_fill_brewer(palette = "Dark2") +
        labs(x = "order", y = "total count")
        P1 <- P1 + geom_bar(stat = "identity")  + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 0))
        P1
        dev.off()

  
      
  

      
#body size of arthropods
      #pull p-values from ANOVAs of body size within order, which are less than 0.05?
      lmp <- function (modelobject) {
        if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
        f <- summary(modelobject)$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        attributes(p) <- NULL
        return(p)
      }
      
      keepers <- list()
      for(i in unique(insectD$order)){ 
        temp <- insectD[insectD$order == i , ]
        temp2 <- lm(length.mm ~ treatment, data = temp)
        keepers[i] <- lmp(temp2)
        }
      which(keepers<0.05)
      
      #do the same with a ks.test for L and U
      
      kTest <- ks.test(
        insectD[insectD$order == "Aranea" & insectD$treatment == "L" , "length.mm"] , 
        insectD[insectD$order == "Aranea" & insectD$treatment == "U" , "length.mm"] )

      keepersKS <- list()
      for(i in unique(insectD$order)){ 
      kTest <- ks.test(
      insectD[insectD$order == i & insectD$treatment == "L" , "length.mm"] , 
      insectD[insectD$order == i & insectD$treatment == "U" , "length.mm"] )
      
      keepers$p.value[i] <- kTest$p.value
      keepers$Dstat[i] <- kTest$statistic
      }
      
    keepersKSnames  <-  which(keepers$p.value<0.05)
    tax.order<-keepers[keepersKSnames]
      #plot densities of taxa

    with(insectD , tapply(length.mm , list(order, treatment) , median))[,1] / 
      with(insectD , tapply(length.mm , list(order, treatment) , median))[,2] 
    
    #a function to plot lit and unlit body size for a taxa
      plot.order.cdf <- function(tax.order){
      plot(ecdf(insectD[insectD$order == tax.order & 
                          insectD$treatment == "L" , "length.mm"]) , 
           lwd = 3,cex.lab = 1.5, cex.axis = 1.5,col = "#E69F00",verticals = TRUE,pch = 16, cex = 0.01 , main = paste(tax.order,"body length")
           ,xlab = "body length (mm)", ylab = "ECDF")
      
      lines(ecdf(insectD[insectD$order == tax.order & 
                          insectD$treatment == "U" , "length.mm"]) , 
           lwd = 3,col = "#999999",verticals = TRUE,pch = 16, cex = 0.01 )
      legend("bottomright", lwd = c(3,3) , col = c("#E69F00","#999999"), legend = c("Lit", "Unlit"))
      }
      
      

      for(i in names(keepersKSnames)){
          pdf(file = paste("body size dist ", i, ".pdf", sep = ""), width = 10, height = 6)
          plot.order.cdf(i)
          dev.off()
      }
      
      with(insectD , tapply(length.mm , list(order, treatment), mean ))
      with(insectD , tapply(length.mm , list(order, treatment), se ))

#predation experiment
      
      #day
       M7 <-  glm(cbind( predExpD$flyStartDay - predExpD$flyEndDay, predExpD$flyEndDay) ~ 
              predExpD$Treatment,
            family = "binomial")      
       M70 <-  glm(cbind( predExpD$flyStartDay - predExpD$flyEndDay, predExpD$flyEndDay) ~ 
                    1,
                  family = "binomial")
       summary(M7)
      anova(M7 , M70, test = "Chisq")
      
      #night
      M8 <-  glm(cbind( predExpD$flyStartNight - predExpD$flyEndNight, predExpD$flyEndNight) ~ 
                   predExpD$Treatment,
                 family = "binomial")      
      M80 <-  glm(cbind( predExpD$flyStartNight - predExpD$flyEndNight, predExpD$flyEndNight) ~ 
                    1,
                  family = "binomial")   
      anova(M8 , M80, test = "Chisq")
      summary(M8)
      #plot
      
      predExpD$propNightEaten <-  (predExpD$flyStartNight - predExpD$flyEndNight )/predExpD$flyStartNight
      predExpD$propDayEaten <-  (predExpD$flyStartDay - predExpD$flyEndDay )/predExpD$flyStartDay
      with(predExpD , tapply(propDayEaten , list(Treatment), mean ))
      
      with(predExpD , tapply(propNightEaten , list(Treatment), mean ))[1]/
        with(predExpD , tapply(propNightEaten , list(Treatment), mean ))[2]
      

      P2 <- ggplot(predExpD , aes(x = Treatment, y = propNightEaten))
      P2<-P2 + geom_bar(stat = "summary", fun.y = "mean", fill = c("#E69F00", "#999999"), color = "white") +
        stat_summary(fun.data = mean_se,geom = "errorbar", width = 0.2) +
        theme_dark()+ 
        scale_y_continuous(limits = c(0, 1))+
        labs(x = "treatment", y = "proportion flies eaten") + ggtitle("Night") +theme(aspect.ratio = 1)  + scale_x_discrete(labels = c("ALAN" , "control"))

      
      
      P1 <- ggplot(predExpD , aes(x = Treatment, y = propDayEaten))
      P1 <- P1 + geom_bar(stat = "summary", fun.y = "mean", fill = c( "#E69F00", "#999999")) +
        stat_summary(fun.data = mean_se,geom = "errorbar", width = 0.2) +
        scale_y_continuous(limits = c(0, 1))+
        labs(x = "treatment", y = "proportion flies eaten") + ggtitle("Day")+theme(aspect.ratio = 1) + scale_x_discrete(labels = c("ALAN" , "control")) +
        theme_bw()
      pdf(file = "../figures/Fig 5 - night day predation.pdf", width = 8, height=4)
      grid.arrange(P1 , P2, ncol = 2)
      dev.off()
      
      #####



