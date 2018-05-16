#### A Start

library(survival)
mydata<-surdf
## Step (1)
## Create the Survival Object
## we need status = 0 --> no event and status = 1 --> event happened
recodestatus<-function(x){
  if(x==1){rs=0} ## no event / censored
  if(x==2){rs=1} ## event happened
  return(rs)
}
for(i in 1:length(mydata$status)){
  mydata$recodedStatus[i]<-recodestatus(mydata$status[i])
}
mySurv<-Surv(time=mydata$Survival.months, event = mydata$Survival.status)
class(mySurv)
head(mySurv)


####### B

## single survival curve: no comparisons
myfit<-survfit(mySurv~1) ## signle curve for all patients in 
## the dataset
myfit
median(mydata$Survival.months,na.rm = TRUE)


### Median survival is the time at which the survivorship 
### function equals 0.5.
plot(myfit)
plot(myfit, conf.int = "none")
abline(h=0.5)
abline(v=52.3)


## specify predictor variable in the formula
myfit<-survfit(mySurv~mydata$pred)
myfit
plot(myfit)
table(mydata$pred)
# 1= Male, 2= Female
plot(myfit, col=c("red","blue","green")) ## red = Male, Blue= female
plot(myfit, conf.int = "both", col=c("red","blue"))
plot(myfit, col=c("red","blue"))
plot(myfit, col=c("red","blue","green"), mark=3) ## mark.time=T marked at 
## each censoring time
legend("topright", c("iCluster1","iCluster2",'iCluster3'), col=c("red","blue","green"), lty=1)
abline(h=0.5)
abline(v=23, col="red")
abline(v=60.5, col="blue")
abline(v=51.5, col="green")

## Now we see that survival of females is better, 
## Q: Is it better by chance, or statistically significant?
survdiff(mySurv~mydata$pred)

###  plot the inverse of a survival function

plot(myfit, fun="event", col=c("red","blue","green"), mark=3)


### End B

### C

coxph(mySurv~mydata$Gender+mydata$Age+mydata$t)

#### End C



currentClin = clinER[, c("OS_dmfs", "OSbin_dmfs", "ER_final", "Age")]
row.names(currentClin) = clinER$SampleID


## -- exprData: genes in cols and sample names in rows
## -- survData: sample names in rows, and survival and covariate information in cols
## -- timeCol: column name that has time for survival
## -- eventCol: column name that has event status for survival
## -- column names for teh covariates to be included

survivalTest = function(exprData = as.matrix(uncorExpr[, -6]),
                        survData = currentClin,
                        timeCol = "OS_dmfs",
                        eventCol = "OSbin_dmfs",
                        covarCol = c("Age", "ER_final"),
                        plotKM = F){
  
  ## remove samples without annotation:
  survData = survData[complete.cases(survData), ]
  
  ## subset expression based on those samples that have OS annotation:
  exprData = exprData[row.names(survData), ]
  
  
  ## calculate median values for all genes and add that to the expression data:
  medExpr = apply(exprData, 2, median)
  exprData = rbind(exprData , medExpr)
  
  ## replace expression values using 
  newData = apply(exprData, 2, function(x){
    sapply(x[-length(x)], function(y){
      ifelse(y > x[length(x)], "High", "Low")
    })
  }) 
  
  row.names(newData) = row.names(exprData[- nrow(exprData),])
  newData = data.frame(newData, check.names = F)
  
  
  coxStats = foreach(g = 1:ncol(newData), .packages="survival", .combine = "rbind") %do% {
    
    coxFit = coxph(Surv(survData[,timeCol], survData[,eventCol]) ~ newData[,g] + 
                     survData[,covarCol[1]] + 
                     survData[,covarCol[2]], 
                   data = survData)
    
    coxRes = c(summary(coxFit)$coef[1,], 
               summary(coxFit)$conf.int[1,2:4], 
               summary(coxFit)$concordance, 
               summary(coxFit)$sctest)
    return(coxRes)
    
  }
  row.names(coxStats) = colnames(newData)
  
  return(as.data.frame(coxStats))
}




##------------------------------------------------------- Example:

uncorStats = survivalTest(exprData = mat,
                          survData = currentClin,
                          timeCol = "OS_dmfs",
                          eventCol = "OSbin_dmfs",
                          covarCol = c("Age", "ER_final"))
