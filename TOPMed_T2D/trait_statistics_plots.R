## Phenotype Plot Statistics
## Alisa Manning (amanning@broadinstitute.org)
## 5 August 2016


########## LOAD LIBRARIES
library(ggplot2)
#install.packages("e1071",repos="http://cran.us.r-project.org/")
library(e1071)

########## FUNCTIONS

get.stats <- function(colname, grouping.label, dataset,sex.column="sex",covariates=c()) {
  print(grouping.label)
  col <- dataset[,colname]
  sex <- dataset[,sex.column]
  print(table(sex))
  sex.labels=names(table(sex))
  print(sex.labels)
  f.names <- c("mean","sd","min","max","skewness","kurtosis")
  to.return <- c()
  to.return.cov <- c()
  
  print(table(!is.na(col)))
  print(head(dataset))
  if( sum(!is.na(col))>0) {
    to.return <- c(column=colname,grouping=grouping.label,all.N=length(col),
     				N=format(sum(!is.na(col))), sapply(f.names,function(f.name){f<-get(f.name);return(f.name=f(col,na.rm=T))}),
     				N.sex1=format(sum(!is.na(col[which(sex==sex.labels[1])]))), sapply(f.names,function(f.name){f<-get(f.name);return(f.name=f(col[which(sex==sex.labels[1])],na.rm=T))}),
     				N.sex2=format(sum(!is.na(col[which(sex==sex.labels[2])]))), sapply(f.names,function(f.name){f<-get(f.name);return(f.name=f(col[which(sex==sex.labels[2])],na.rm=T))})
     				)
     for(cov in covariates) {
     		cov.col <- dataset[,cov]
     		cov.col[which(is.na(col))] <- NA
     		
     	    to.return.cov.1 <- c(column=paste(colname,cov,sep="_cov_"),grouping=grouping.label,all.N=length(cov.col),
     				N=format(sum(!is.na(cov.col))), sapply(f.names,function(f.name){f<-get(f.name);return(f.name=f(cov.col,na.rm=T))}),
     				N.sex1=format(sum(!is.na(cov.col[which(sex==sex.labels[1])]))), sapply(f.names,function(f.name){f<-get(f.name);return(f.name=f(cov.col[which(sex==sex.labels[1])],na.rm=T))}),
     				N.sex2=format(sum(!is.na(cov.col[which(sex==sex.labels[2])]))), sapply(f.names,function(f.name){f<-get(f.name);return(f.name=f(cov.col[which(sex==sex.labels[2])],na.rm=T))})
     				)
     				
     		to.return.cov <- rbind(to.return.cov,to.return.cov.1)
     }
     
  } else {
    to.return <- c(column=colname,grouping=grouping.label,all.N=length(col),N=format(sum(!is.na(col))),rep(NA,6),
    N.sex1=format(sum(!is.na(col[which(sex==sex.labels[1])]))),rep(NA,6),
    N.sex2=format(sum(!is.na(col[which(sex==sex.labels[2])]))),rep(NA,6))
  } 
  names(to.return) <- c("colname","Group","N.total","N.nonmissing",f.names,paste(c("N.nonmissing",f.names),paste(".sex",sex.labels[1],sep=""),sep=""), paste(c("N.nonmissing",f.names),paste(".sex",sex.labels[2],sep=""),sep=""))
  if(!is.null(to.return.cov)) {
  	to.return <- rbind(to.return,to.return.cov)
  }
  return(to.return)
}




##################################################

############### 
# Generate sample statistics and plots for each QT




stats <- c()

for(group in grouping) {

for(trait in glycemic.traits) {
	trait.covs <- c(age.map[[trait]],bmi.map[[trait]])
	
	if(trait==hba1c.trait.name) {
		trait.covs <- c(trait.covs,hba1c.covs.map[["HbA1c.FG"]], hba1c.covs.map[["HbA1c.2hr"]])
	}
	

	for(level in unique(in.ped[,group])) {
		stats <- rbind(stats,
					get.stats(trait,paste("Group",group,"Level",level,sep="_"),in.ped[which(in.ped[,group]==level & in.ped[,t2d.map[[trait]]]==0),],covariates=trait.covs),						
					get.stats(trait,paste("Group",group,"Level",level,"SEQUENCED",sep="_"),in.ped[which(in.ped[,"sequenced"]==1 & in.ped[,group]==level & in.ped[,t2d.map[[trait]]]==0),],covariates=trait.covs),
					get.stats(trait,paste("Group",group,"Level",level,"SEQUENCED_FRAM_WGS",sep="_"),in.ped[which(in.ped[,"sequenced.cohort"]=="FRAM_WGS" & in.ped[,group]==level & in.ped[,t2d.map[[trait]]]==0),],covariates=trait.covs),
					get.stats(trait,paste("Group",group,"Level",level,"SEQUENCED_EOAF_WGS",sep="_"),in.ped[which(in.ped[,"sequenced.cohort"]=="EOAF_WGS" & in.ped[,group]==level & in.ped[,t2d.map[[trait]]]==0),],covariates=trait.covs)
				)
	}
}
}
	###
write.table(stats,paste(cohort,".glycemic.distributions.",paste(strsplit(date(),split=" +")[[1]][c(3,2,5)],collapse=""),".csv",sep=""),row.names=F,quote=F,sep=",")


pdf(paste(cohort,".glycemic.distributions.",paste(strsplit(date(),split=" +")[[1]][c(3,2,5)],collapse=""),".pdf",sep=""),width=9,height=6)

for(trait in glycemic.traits) {
    print(ggplot(in.ped[which(in.ped[,t2d.map[[trait]]]==0 & !is.na(in.ped[,trait])),],aes_string(x=trait)) 
    			+ geom_density() 
    			+ geom_density(aes(col=factor(sequenced.cohort))) 
    			+ facet_grid(sex ~ cohort, labeller="label_both") 
    			+ theme_classic())
    
      
        print(ggplot(in.ped[which(in.ped[,t2d.map[[trait]]]==0 & !is.na(in.ped[,trait])),],aes_string(x=trait)) 
    			+ geom_histogram() 
    			+ facet_grid(sex ~ cohort, labeller="label_both") 
    			+ theme_classic())
  
    
        print(ggplot(in.ped[which(in.ped[,t2d.map[[trait]]]==0 & !is.na(in.ped[,trait])),],aes_string(x=trait)) 
    			+ geom_histogram(aes(col=factor(sequenced.cohort))) 
    			+ facet_grid(sex ~ cohort, labeller="label_both") 
    			+ theme_classic())
    			
    print(ggplot(in.ped[which(in.ped[,t2d.map[[trait]]]==0 & !is.na(in.ped[,trait])),],aes_string(y=trait)) 
    			+ geom_boxplot(aes(x=factor(sequenced.cohort),col=factor(sequenced.cohort))) 
    			+ facet_grid(sex ~ cohort, labeller="label_both") 
    			+ theme_classic())
    		


}
dev.off()



############################################
##### T2D distributions


pdf(paste(cohort,".t2d.distributions.",paste(strsplit(date(),split=" +")[[1]][c(3,2,5)],collapse=""),".pdf",sep=""),width=11,height=5.5)

## note: changed from in.t2d[which(in.t2d[,"T2D"] == level),] to in.t2d[which(in.t2d[,"T2D"] %in% level),] because there are some 
## individuals in ped file where T2D values are NA
#
#> table(in.t2d$T2D,useNA="always")
#
#    0     1  <NA>
#12415  1905    81

t2d.stats <- c()
for(trait in t2d.continuous.covariates) {
		print(trait)
			t2d.stats <- rbind(t2d.stats,
						get.stats(trait,"All",in.t2d),			
						get.stats(trait,"SEQUENCED",in.t2d[which(in.t2d[,"sequenced"]==1),]))
		
		for(level in unique(in.ped[,"T2D"])) {
			print(paste("T2D level",level))
			t2d.stats <- rbind(t2d.stats,
						get.stats(trait,paste("Group","T2D","Level",level,sep="_"),in.t2d[which(in.t2d[,"T2D"] %in% level),]),						
						get.stats(trait,paste("Group","T2D","Level",level,"SEQUENCED",sep="_"),in.t2d[which(in.t2d[,"sequenced"]==1 & in.t2d[,"T2D"] %in% level),])
						)
		}
		### cohort
		for(cohort.1 in unique(in.t2d$cohort)) {
			t2d.stats <- rbind(t2d.stats,
						get.stats(trait,paste("Cohort",cohort.1,sep="_"),in.t2d[which(in.t2d[,"cohort"] == cohort.1),]),			
						get.stats(trait,paste("Cohort",cohort.1, "SEQUENCED",sep="_"),in.t2d[which(in.t2d[,"cohort"] == cohort.1 & in.t2d[,"sequenced"]==1),]))
		
			for(level in unique(in.ped[,"T2D"])) {
				print(paste("T2D level",level))
				t2d.stats <- rbind(t2d.stats,
						get.stats(trait,paste("Cohort",cohort.1,"Group","T2D","Level",level,sep="_"),in.t2d[which(in.t2d[,"cohort"] == cohort.1 & in.t2d[,"T2D"] %in% level),]),						
						get.stats(trait,paste("Cohort",cohort.1,"Group","T2D","Level",level,"SEQUENCED",sep="_"),in.t2d[which(in.t2d[,"cohort"] == cohort.1 & in.t2d[,"sequenced"]==1 & in.t2d[,"T2D"] %in% level),])
						)
		}
						
	}

   print(ggplot(in.t2d[which(!is.na(in.t2d[,trait])),],aes_string(x=trait)) 
    			+ geom_density() 
    			+ geom_density(aes(col=factor(sequenced.cohort))) 
    			+ facet_grid(sex * T2D ~ cohort, labeller="label_both") 
    			+ theme_classic())
    

   print(ggplot(in.t2d[which(!is.na(in.t2d[,trait])),],aes_string(x=trait)) 
    			+ geom_histogram() 
    			+ facet_grid(sex * T2D ~ cohort, labeller="label_both") 
    			+ theme_classic())
 
    print(ggplot(in.t2d[which(!is.na(in.t2d[,trait])),],aes_string(x=trait)) 
    			+ geom_histogram(aes(col=factor(sequenced.cohort))) 
    			+ facet_grid(sex * T2D ~ cohort, labeller="label_both") 
    			+ theme_classic())

  print(ggplot(in.t2d[which(!is.na(in.t2d[,trait])),],aes_string(y=trait)) 
    			+ geom_boxplot(aes(x=factor(sequenced.cohort),col=factor(sequenced.cohort))) 
    			+ facet_grid(sex * T2D ~ cohort, labeller="label_both") 
    			+ theme_classic())
}

dev.off()

write.table(t2d.stats,paste(cohort,".t2d.distributions.",paste(strsplit(date(),split=" +")[[1]][c(3,2,5)],collapse=""),".csv",sep=""),row.names=F,quote=F,sep=",")


