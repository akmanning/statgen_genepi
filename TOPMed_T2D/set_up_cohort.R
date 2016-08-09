
#### INPUT DATA

############# Glycemic trait ped file

# label for cohort
cohort <- "FHS"

in.ped <- read.table("TOPMed.glycemic.traits.8Aug2016.ped",header=T,as.is=T)
in.ped$logFI <- log(in.ped$FastingInsulin)

#in.ped$cohort.T2D <- paste(in.ped$cohort,in.ped$T2D,sep=".")



# Change column names to match the input file
glycemic.traits <- c("FastingGlucose","FastingInsulin","logFI","HbA1c")

# the column name used for HbA1c in this data set
hba1c.trait.name <- "HbA1c"

# specify any grouping variables.

in.ped$cohort.sequenced.cohort <- paste(in.ped$cohort,in.ped$sequenced.cohort,sep=".")

grouping <- c("cohort") #"cohort.sequenced.cohort") # ,"sequenced.cohort" ,"T2D","cohort.T2D")


# provide label for T2D status for each glycemic trait
t2d.map <- list("FastingGlucose"="T2D_FG", "FastingInsulin"="T2D_FI","logFI"="T2D_FI", "HbA1c"="T2D_HbA1c")
age.map <- list("FastingGlucose"="age_FG", "FastingInsulin"="age_FI","logFI"="age_FI", "HbA1c"="age_HbA1c")
bmi.map <- list("FastingGlucose"="BMI_FG", "FastingInsulin"="BMI_FI","logFI"="BMI_FI", "HbA1c"="BMI_HbA1c")
hba1c.covs.map <- list("HbA1c.FG"="FastingGlucose_HbA1c","HbA1c.2hr"="Gluc2Hr_HbA1c")

glycemic.traits %in% colnames(in.ped) ## all list items should be TRUE
unlist(t2d.map) %in% colnames(in.ped)
unlist(age.map) %in% colnames(in.ped)
unlist(bmi.map) %in% colnames(in.ped)
unlist(hba1c.covs.map) %in% colnames(in.ped)

############ T2D ped file

in.t2d <- read.csv("/restricted/projectnb/glycemic/phenotypes/T2D_G1ex27_G2ex9_g3ex2_noOGTT_05092016.csv",header=T,as.is=T)
t2d.continuous.covariates <- c("last_exam_age","last_exam_BMI","last_exam_FG","T2D_age","T2D_BMI")
t2d.categorical.covariates <- c("sex","last_exam_T2D_treatment")

## remove people with T2D==NA ? 
#in.t2d <- in.t2d[which(!is.na(in.t2d[,"T2D"])),]

in.seq <- read.csv("/restricted/projectnb/glycemic/ctliu/TOPMed/TOPMED_ID.csv",header=T,as.is=T)
in.seq$sequenced <- ifelse(in.seq$EOAF_WGS==1 | in.seq$FRAM_WGS==1, 1,0)
in.seq$sequenced.cohort <- apply(in.seq,1,function(x){if(x["EOAF_WGS"]==1) {return("EOAF_WGS")} else if (x["FRAM_WGS"]==1) {return("FRAM_WGS")} else {return("Not_sequenced")}})

table(in.t2d$id %in% in.seq$framid)
dim(in.t2d)
dim(in.seq)
in.t2d <- merge(in.t2d,in.seq,all.x=T,by.x="id",by.y="framid")

in.t2d$cohort <- factor(in.t2d[,"idtype"], levels=c(0,1,3),labels=c("Original","Offspring","Gen3"))

#source("trait_statistic_plot.R")

