# The R code is written in collaborative work of Dr. Ji Eun Park and Prof. Seo Young Park 
# please contact jieunp@gmail.com if you have further question.
# These are codes for radiogenomics analysis "prediction of core signaling pathway in IDH-wildtype glioblastoma"

## Line 5-129: Feature selection via Student's t-test with false discovery rate correction
## Line 134-233: Feature selection via LASSO penalization and calculate AUC for each genetic mutation
## Line 234-329: Feature selection via Random Forest and find top 5 important features.
## Line 330- : Calculate diagnostic performance

#### First step is univariate analysis 

rm(list=ls())
setwd("C:/Radiogenomics/") #write down where the feature exists


rm(list=ls())


require(plyr)
require(dplyr)

#### Read the feature data in
t1<-read.csv("Feature_T1_Z.csv",header = FALSE)
flair<-read.csv("Feature_FLAIR_Z.csv",header = FALSE)
adc<-read.csv("Feature_ADC.csv",header = FALSE)
cbv<-read.csv("Feature_DSC.csv",header = FALSE)

t1_t<-read.csv("Test_Feature_T1_Z.csv",header = FALSE)
flair_t<-read.csv("Test_Feature_FLAIR_Z.csv",header = FALSE)
adc_t<-read.csv("Test_Feature_ADC.csv",header = FALSE)
cbv_t<-read.csv("Test_Feature_DSC.csv",header = FALSE)

t1.n<-cbind(t1,t1_t)
flair.n<-cbind(flair,flair_t)
adc.n<-cbind(adc,adc_t)
cbv.n<-cbind(cbv,cbv_t)

all<-rbind(t1.n,flair.n,adc.n,cbv.n)
xall<-t(all)

colnames(xall)<-c(paste("t1.", 1:nrow(t1.n), sep=""), 
                  paste("flair.", 1:nrow(flair.n), sep=""),
                  paste("adc.", 1:nrow(adc.n), sep=""),
                  paste("cbv.", 1:nrow(cbv.n), sep="")  )

#summary(xall)
#### Read the genomic data in
g.data1<-read.csv("genomic_train.csv",header=TRUE)
g.data2<-read.csv("genomic_test.csv",header=TRUE)

g.data<-rbind(g.data1,g.data2)
#### Read the genomic data in

summary(g.data)

# number of genes 
no.g<-16 #IDH, 1p19q include 16,17

# Student's t-test with false discovery rate to select features 

selected.f<-as.list(rep(NA, no.g))
names(selected.f)<-colnames(g.data)[2:17] #IDH, 16, 1p19q 17
# number of selected features
no.s.f<-rep(NA, no.g)

# Instead, we decided to use p-value to screen the features
for(i in 1:no.g){
  temp<-p.adjust(apply(xall, 2, t.fun, gg=g.data[, i+1]), method="fdr")
  selected.f[[i]]<-which(temp<0.5)
  no.s.f[i]<-length(selected.f[[i]])
}
no.s.f


egfr<-c(selected.f[["EGFR"]])
pdgf<-c(selected.f[["PDGFRA"]])
pi3k<-c(selected.f[["PI3K"]])
pten<-c(selected.f[["PTEN"]])
nf1<-c(selected.f[["NF1"]])
mdm2<-c(selected.f[["MDM2"]])
p53<-c(selected.f[["P53"]])
cdk4<-c(selected.f[["CDK4"]])
cdkn<-c(selected.f[["CDKN2AB"]])
rb<-c(selected.f[["Rb1"]])
ccnd<-c(selected.f[["CCND2"]])
IDH<-selected.f[["IDHm"]]

length(egfr)
length(pdgf)
length(pi3k)
length(pten)
length(nf1)
length(mdm2)
length(p53)
length(cdk4)
length(cdkn)
length(rb)
length(ccnd)
length(IDH)


### data of the selected features
#### center and scale?
#x.data.rtk<-scale(xall[, RTK.related])
x.data.egfr<-xall[,egfr]
x.data.pdgf<-xall[,pdgf]
x.data.pi3k<-xall[,pi3k]
x.data.pten<-xall[,pten]
x.data.nf1<-xall[,nf1]
x.data.mdm2<-xall[,mdm2]
x.data.p53<-xall[,p53]
x.data.cdk4<--xall[,cdk4]
x.data.cdkn<-xall[,cdkn]
x.data.rb<-xall[,rb]
x.data.ccnd<-xall[,ccnd]
x.data.idh<-xall[,IDH]

setwd("C:/Radiogenomics/result_t_test")


write.table(x.data.egfr, file="egfr_selected.csv",sep=",",row.names=F, col.names=TRUE)
write.table(x.data.pdgf, file="pdfg_selected.csv",sep=",",row.names=F, col.names=TRUE)
write.table(x.data.pi3k, file="pi3k_selected.csv",sep=",",row.names=F, col.names=TRUE)
write.table(x.data.pten, file="pten_selected.csv",sep=",",row.names=F, col.names=TRUE)
write.table(x.data.nf1, file="nf1_selected.csv",sep=",",row.names=F, col.names=TRUE)
write.table(x.data.mdm2, file="mdm2_selected.csv",sep=",",row.names=F, col.names=TRUE)
write.table(x.data.p53, file="p53_selected.csv",sep=",",row.names=F, col.names=TRUE)
write.table(x.data.cdk4, file="cdk4_selected.csv",sep=",",row.names=F, col.names=TRUE)
write.table(x.data.cdkn, file="cdkn_selected.csv",sep=",",row.names=F, col.names=TRUE)
write.table(x.data.rb, file="rb_selected.csv",sep=",",row.names=F, col.names=TRUE)
write.table(x.data.ccnd, file="ccnd_selected.csv",sep=",",row.names=F, col.names=TRUE)
write.table(x.data.idh,file="idh_selected.csv",sep=",",row.names=F, col.names=TRUE) 



####### The second step is LASSO based feature selection for each gene
# change xdat at each time by removing "#"

xdat<-read.csv("egfr_selected.csv")
#xdat<-read.csv("pdfg_selected.csv")
#xdat<-read.csv("pi3k_selected.csv")
#xdat<-read.csv("pten_selected.csv")
#xdat<-read.csv("nf1_selected.csv") 
#xdat<-read.csv("mdm2_selected.csv")
#xdat<-read.csv("p53_selected.csv") 
#xdat<-read.csv("cdk4_selected.csv") 
#xdat<-read.csv("cdkn_selected_n.csv") 
#xdat<-read.csv("rb_selected.csv")  
#xdat<-read.csv("idh_selected.csv")

dim(xdat)



#### Read the genomic data in

g.data<-read.csv("genomic_train_and_test_order_120.csv",header=T)

summary(g.data)

y<-g.data$EGFR
#y<-g.data$PDGFRA
#y<-g.data$PI3K
#y<-g.data$PTEN
#y<-g.data$NF1
#y<-g.data$MDM2
#y<-g.data$P53
#y<-g.data$CDK4
#y<-g.data$CDKN2AB
#y<-g.data$Rb1
#y<-g.data$IDHm

y<-as.factor(y)

require(dplyr)


## Additional analysis) xdat was further separated as anatomic, diffusion, and perfusion imaging during the revision
#x.ana<-xdat %>% select(starts_with("t1."), starts_with("flair."))
#x.ana<-xdat %>% select(starts_with("adc."))
#x.ana<-xdat %>% select(starts_with("cbv."))


x.train <-data.matrix(xdat[1:85,]) # change xdat as x.ana for additional analysis
x.test <- data.matrix(xdat[86:120,])

yy.train <-data.matrix(y[1:85])
yy.test <-data.matrix(y[86:120])


#LASSO based selection and diagnostic performance using 3-fold cross validation 

set.seed(600) 
require(glmnet)
require(pROC)

glmnet.obj <- cv.glmnet(x.train,yy.train, family=c("binomial"),nfolds=3)
train.p<-predict(glmnet.obj,newx=x.train,s="lambda.min",type="response")
auc.tr<-roc(yy.train~as.vector(train.p))  
test.p<-predict(glmnet.obj, 
                newx=x.test, s="lambda.min", type="response")
auc.tst<-roc(yy.test~as.vector(test.p)) 

auc.tr
ci(auc.tr)
coords(auc.tr, "best", ret=c("threshold", "sens", "spec", "ppv", "npv", "accuracy"))

auc.tst
ci(auc.tst)
coords(auc.tst, "best", ret=c("threshold", "sens", "spec", "ppv", "npv", "accuracy"))


coef.mat<-matrix(NA, 1+ncol(x.ana), 1)
coef.mat[1:nrow(coef.mat),1]<-as.vector(coef(glmnet.obj, s="lambda.min"))       
rownames(coef.mat)<-c(("intercept"),colnames(x.ana))

selected.mat<-coef.mat[(coef.mat!=0)]
names(selected.mat)<-rownames(coef.mat)[coef.mat!=0]

xx <-x.ana[,names(selected.mat)[-1]]

# coef. name needs
setwd("C:/Radiogenomics/lasso_results")

write.table(xx,file="coeff_egfr.csv",sep=",")
#write.table(xx,file="coeff_pdfg.csv",sep=",")
#write.table(xx,file="coeff_pi3k.csv",sep=",")
#write.table(xx,file="coeff_pten.csv",sep=",")  
#write.table(xx,file="coeff_nf1.csv",sep=",")
#write.table(xx,file="coeff_mdm2.csv",sep=",")
#write.table(xx ,file="coeff_p53.csv",sep=",")
#write.table(xx,file="coeff_cdk4.csv",sep=",")
#write.table(xx,file="coeff_cdkn.csv",sep=",")
#write.table(xx,file="coeff_rb.csv",sep=",")
#write.table(xx,file="coeff_ccnd.csv",sep=",")
#write.table(xx,file="coeff_idh.csv",sep=",")


#### Third Step: Random Forest Classifier for Feature Selection
##Before this step , you need to make a combined feature file for each pathway

rm(list=ls())

setwd("C:/Radiogenomics/lasso_results")

egfr<- read.csv("coeff_egfr.csv",header=T) 
pdgf<- read.csv("coeff_pdfg.csv",header=T) 
pi3k<- read.csv("coeff_pi3k.csv",header=T) 
pten<- read.csv("coeff_pten.csv",header=T) 
mdm2<- read.csv("coeff_mdm2.csv",header=T) 
p53<- read.csv("coeff_p53.csv",header=T) 
cdk4<- read.csv("coeff_cdk4.csv",header=T) 
cdkn<- read.csv("coeff_cdkn.csv",header=T) 
rb<- read.csv("coeff_rb.csv",header=T) 
ccnd<- read.csv("coeff_ccnd.csv",header=T)

path_rtk <-cbind (egfr,pdfg,pi3k,pten)
path_p53 <-cbind (mdm2,p53)
path_rb <-cbind (cdk4,cdkn,rb,ccnd)

# check and remove the overlapped features in the path_rtk, path_p53, and path_rb then save them

setwd("C:/Radiogenomics/pathway")


write.table(path_rtk,file="coeff_path_RTK.csv",sep=",")
write.table(path_p53,file="coeff_path_P53.csv",sep=",")
write.table(path_rb,file="coeff_path_Rb.csv",sep=",")

##
rm(list=ls())

setwd("C:/Radiogenomics/pathway")

# change the xdat at each time
xdat <- read.csv("coeff_path_RTK.csv",header=T) 
#xdat <- read.csv("coeff_path_P53.csv",header=T) 
#xdat <- read.csv("coeff_path_Rb.csv",header=T)


summary(xdat)

setwd("C:/Radiogenomics/")
g.data<-read.csv("genomic_train_and_test_order_120.csv",header=T)

summary(g.data)

y<-g.data$pathway_RTK
#y<-g.data$pathway_p53
#y<-g.data$pathway_Rb
#y<-g.data$IDHm

xtr <-data.matrix(xdat[1:85,])
xtst <- data.matrix(xdat[86:120,])
ytr<-data.matrix(y[1:85])
ytst<-as.factor(y[86:120])


#########random Forest

require(randomForest)

set.seed(100)

#train data 
train.data <-data.frame(ytr,xtr)
train.data$ytr<-as.factor(train.data$ytr)

#test data
test.data <-data.frame(ytst,xtst)
test.data$ytst<-as.factor(test.data$ytst)

fitrf <-randomForest(ytr ~.,train.data, importance=TRUE, ntree =100)

plot(fitrf)

print(fitrf)

summary(fitrf)

round(importance(fitrf),2)

varImpPlot(fitrf)

importance(fitrf)


measure<-data.frame(importance(fitrf))
measure$Vars<-row.names(measure)
arrange(measure,desc(measure$MeanDecreaseGini))
arrange(measure,desc(measure$MeanDecreaseAccuracy))

#Diagnostic performance 
# RTK
fit<-glm(ytr~t1.xx+t1.xx+flair.xx+adc.xx+cbv.xx, family = binomial(link = "logit"),data=train.data)  # this is an example. Write down the selected top 5 features from the Random Forest
fit.tst<-glm(ytst~t1.xx+t1.xx+flair.xx+adc.xx+cbv.xx, family = binomial(link = "logit"),data=test.data) # this is an example.


#P53
fit<-glm(ytr~t1.930+t1.412+t1.439+t1.540+t1.951, family = binomial(link = "logit"),data=train.data) # this is an example.
fit.tst<-glm(ytst~t1.930+t1.412+t1.439+t1.540+t1.951, family = binomial(link = "logit"),data=test.data) # this is an example.

#Rb
fit<-glm(ytr~t1.743+t1.1443+t1.307+t1.1427+flair.910, family = binomial(link = "logit"),data=train.data)# this is an example.
fit.tst<-glm(ytst~t1.743+t1.1443+t1.307+t1.1427+flair.910, family = binomial(link = "logit"),data=test.data)# this is an example.


#Diagnostic performance 

pred=predict(fit)
roc1=roc (ytr ~ pred)

plot (roc1)

roc1
ci(roc1)
coords(roc1, "best", ret=c("threshold", "sens", "spec", "ppv", "npv", "accuracy"))
cutoff.value<-coords(roc1, "best", ret=c("threshold", "sens", "spec", "ppv", "npv", "accuracy"))[["threshold"]]
predicted.outcome<-as.factor(ifelse(pred>=cutoff.value, 1, 0))
confusionMatrix(table(predicted.outcome, ytr))

pred.tst=predict(fit.tst)
roc2=roc (ytst ~ pred.tst)

plot (roc2)

roc2
ci(roc2)
coords(roc2, "best", ret=c("threshold", "sens", "spec", "ppv", "npv", "accuracy"))
cutoff.value.tst<-coords(roc2, "best", ret=c("threshold", "sens", "spec", "ppv", "npv", "accuracy"))[["threshold"]]
predicted.outcome.tst<-as.factor(ifelse(pred>=cutoff.value.tst, 1, 0))
confusionMatrix(table(predicted.outcome.tst, ytst))

