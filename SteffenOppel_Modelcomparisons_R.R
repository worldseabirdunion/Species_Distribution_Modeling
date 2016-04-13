##############################################################################################
# Electronic Supplement to "Comparison of five modelling techniques to			     #
# predict the spatial distribution and abundance of seabirds"				     #
#											     #
# S. Oppel, A. Meirinho, I. Ramirez, B. Gardner, A. O’Connell, P. Miller, M. Louzao          #
##############################################################################################

# R Program code written by Steffen Oppel, 15 October 2011, with assistance by Maite Louzao
# for questions and comments please contact: steffen.oppel@rspb.org.uk
# Models run in R 2.13.1 on Windows 7 (x64 bit) and in R 2.11.1 on Windows2000 PC



#########################  LOADING REQUIRED FUNCTIONS AND PACKAGES ###########################
#########################  in addition: need Maxent v. 3.3.3 installed from:    ##############
#########################  http://www.cs.princeton.edu/~schapire/maxent/        ##############

setwd("C:/whatever/MyModels")                            ### ADJUST TO PERSONAL DIRECTORY ####
library(PresenceAbsence)
library(gbm)
library(foreign)
library(Design)
library(rJava)
library(randomForest)
library(MuMIn)
library(ROCR)
library(RODBC)						 ### NOTE THAT THIS WILL NOT WORK IN A (x64-bit system)
library(spdep)
library(mgcv)
library(lmodel2)
library(Hmisc)
library(dismo)						 
source("plots.R")					 ### custom functions distributed by Phillips, S. J., and J. Elith. 2010. POC plots: calibrating species distribution models with presence-only data. Ecology 91:2476-2484
source("brt.functions.R")                                ### custom functions distributed by Elith, J., J. R. Leathwick, and T. Hastie. 2008. A working guide to boosted regression trees. Journal of Animal Ecology 77:802-813






#########################  LOADING DATA FROM DATABASE ########################################
########### Please contact the authors if you wish to use the original data  #################
# The column labelled 'SPECIES' contains the presence response variable (0 when absent, 1 when present)
# The column labelled 'COUNT' contains the abundance response variable (i.e. number of birds seen per square km)
# The remaining columns contain predictor variables (SST, CHL-A, etc.)
# If you are running R under a x64-bit system, import data via the 'read.table' command - RODBC not yet supported for 64-bit


channel1<-odbcConnectAccess('SPEA_seabird_data.mdb')
SPEA <- sqlQuery(channel1, "SELECT * FROM SPECIES_all_obs")
odbcClose(channel1)
SPEA$Year<-as.factor(SPEA$Year)
SPEA$season<-as.factor(SPEA$season)




#########################  SETTING TRAINING AND TEST DATA ####################################

t1<-subset(subset(SPEA, LAT>39.7), COUNT>-9999)
t2<-subset(subset(SPEA, LAT<38), COUNT>-9999)
TRAIN<-merge(t1, t2, all=T)
TEST<-subset(subset(subset(SPEA, LAT>=38), LAT<=39.7), COUNT>-9999)
TRAIN$SPECIES<-ifelse(TRAIN$COUNT>0,1,0)
TEST$SPECIES<-ifelse(TEST$COUNT>0,1,0)
remove(t1, t2)



#########################  CORRELATION OF ENVIRONMENTAL VARIABLES ############################

cor.matrix<-as.matrix(TRAIN[,8:18])
cor.matrix<-rcorr(cor.matrix, type="pearson")
P.val<-p.adjust(cor.matrix$P, method = 'bonferroni')
write.table(cor.matrix$r,"Correlation_matrix.csv", sep=",", row.names=F, quote=F)
write.table(P.val,"Correlation_matrix_p_values.csv", sep=",", row.names=F, quote=F)




#########################  TEST FOR AUTOCORRELATION USING MORAN's I AND GEARY's C##########################
################ taken from http://www.ats.ucla.edu/stat/r/faq/morans_i.htm   ################

MoranI<-data.frame(season=c(1:13))
MoranI$statistic<-0
MoranI$p<-0
for (i in 1:13){								### loop to calculate Morans I for each season separately
  b<-subset(TRAIN, Period==i)
  b<-subset(b, COUNT>0)
  SPECIES.coords <- as.matrix(cbind(b$LONG, b$LAT))
  SPECIES.dists<-knn2nb(knearneigh(SPECIES.coords, k= 50, longlat=TRUE))		### adjust k= to the number of nearest neighbours
  coord.list <-make.sym.nb(SPECIES.dists) 
  coord.list <- nb2listw(coord.list,glist=NULL,style="W",zero.policy=FALSE) 
  out<-moran.test(b$COUNT, coord.list, randomisation=TRUE, zero.policy=TRUE, na.action=na.omit, spChk=NULL, adjust.n=TRUE)
  out2<-geary.test(b$BASH, coord.list, randomisation=TRUE, zero.policy=TRUE, spChk=NULL, adjust.n=TRUE)
  MoranI[i,2]<-as.numeric(out$estimate[1])
  MoranI[i,3]<-out$p
  MoranI[i,4]<-as.numeric(out2$estimate[1])
  MoranI[i,5]<-out2$p
}
MoranI
write.table(MoranI, "Spatial_autocorrelation.csv", sep=",", row.names=F)





##############################################################################################
#########################  CONSTRUCTING THE DISTRIBUTION MODELS ##############################
###################### WITH SEASON AND YEAR AS FACTOR PREDICTORS   ###########################
##############################################################################################



#########################  BOOSTED REGRESSION TREES ##########################################

BRT <- gbm.step(data=TRAIN, gbm.x = 3:20, gbm.y = 21, family = "bernoulli", tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.64)




#########################  RANDOM FOREST #####################################################

RF<-randomForest(SPECIES~Year+season+LONG+LAT+BATIMETRY+front_distance+front_dens+gradient+DIST+SST+SSTgrad+SSTchange+SSTanomaly+CHL+CHLgrad+CHLchange+SSH, data=TRAIN, importance=T, replace=T, mtry=2, ntree=1500, na.action=na.omit)




#########################  GLM  ##############################################################

global<-glm(SPECIES~Year+season+LONG+LAT+BATIMETRY+front_distance+front_dens+DIST+I(DIST^2)+SST+SSTgrad+SSTchange+SSTanomaly+CHL+CHLgrad+CHLchange+SSH, data=TRAIN,  family=binomial)
GLM<-stepAIC(global, trace = 0, na.action=na.omit)
#GLM<-dredge(global, beta=F, eval=T, m.max=10,rank="AICc", marg.ex = NULL, trace=F)			### very memory intensive, and yielded same model in our example, not recommended!




#########################  GAM  ##############################################################

detach("package:gam",unload=TRUE)   									### the 'gam' function exists in two packages, and package 'gam' is loaded automatically as dependent of another package. Detach it to call the 'gam' function in package 'mgcv'
GAM<-gam(SPECIES~Year+season+s(LONG,k=5)+s(LAT,k=5)+s(BATIMETRY, k=5)+s(front_distance,k=5)+s(front_dens,k=5)+s(gradient,k=5)+s(DIST, k=5)+s(SST, k=5)+s(SSH, k=5)+s(SSTgrad, k=5)+s(SSTchange, k=5)+s(SSTanomaly, k=5)+s(CHL, k=5)+s(CHLgrad, k=5)+s(CHLchange, k=5), data=TRAIN,  family=binomial, select=TRUE, method='GCV.Cp')



#########################  MAXIMUM ENTROPY  ##################################################
#########  WILL ONLY RUN IF THE 'maxent.jar' file is in the appropriate folder   #############
### NOTE: if maxent throws an error message regarding 'rjava', try the following (with the appropriate file paths for your system):
#.jaddClassPath("C:/Program Files/Java/jre6/lib")
#.jaddClassPath("C:/Program Files/Java/jre6/bin")
#.jaddClassPath("C:/Program Files/R/R-2.13.1/library/dismo/java/maxent.jar")
#.jinit()
# if this does not work, you may have to set the environment variables for Java:
# Desktop-> right click on 'My Computer'->Properties->Advanced->Environment Variables->NEW-> 'Variable name': JAVA_HOME, 'Variable value': C:\Program Files (x86)\Java\jre6


p<-as.vector(TRAIN$SPECIES)
occ<-TRAIN
occ$SPECIES<-NULL
occ$COUNT<-NULL
occ$effort<-NULL
MAX <- maxent(occ, p, args="defaultprevalence=0.1")







#########################  PREDICTION TO TRAINING DATA SET ####################################

TRAIN$BRT_pred<-predict.gbm(BRT, TRAIN, n.trees=BRT$gbm.call$best.trees, type="response")
TRAIN$RF_pred<-predict(RF, TRAIN)
TRAIN$GLM_pred<-predict(GLM, type="response", newdata=TRAIN)
TRAIN$GAM_pred<-as.numeric(predict(GAM, type="response", newdata=TRAIN))
test2<-TRAIN
test2$SPECIES<-NULL
test2$Year<-as.numeric(test2$Year)				### the 'predict' function for maxent requires all numerical variables
test2$season<-as.numeric(test2$season)				### the 'predict' function for maxent requires all numerical variables
test2$COUNT<-NULL
TRAIN$MAX_pred<- predict(MAX, test2, factors=c('season', 'Year'))
rm(test2)

##### ENSEMBLE PREDICTION requires AUC values to calculate the weights for the ensemble prediction
DATA2<-data.frame(TRAIN$LAT, TRAIN$SPECIES, TRAIN$BRT_pred, TRAIN$RF_pred, TRAIN$GLM_pred, TRAIN$GAM_pred, TRAIN$MAX_pred)
Model_evaluation_table_train<-presence.absence.accuracy(DATA2)
weights_t<-(Model_evaluation_table_train[,7]-0.5)/(sum(Model_evaluation_table_train[,7]-0.5))
TRAIN$Ensemble_pred<-((TRAIN$BRT_pred*weights_t[1])+(TRAIN$RF_pred*weights_t[2])+ (TRAIN$GLM_pred*weights_t[3]) + (TRAIN$GAM_pred*weights_t[4]) + (TRAIN$MAX_pred*weights_t[5]))




##############  MODEL VALIDATION AND CALIBRATION OF TRAINING DATA FIT ##########################

Model_calibration_table_train<-data.frame()
DATA2<-data.frame(TRAIN$LAT, TRAIN$SPECIES, TRAIN$BRT_pred, TRAIN$RF_pred, TRAIN$GLM_pred, TRAIN$GAM_pred, TRAIN$MAX_pred, TRAIN$Ensemble_pred)
Model_calibration_table_train<-data.frame("Model"=c("BRT", "RF", "GLM", "GAM", "MAX", "Ensemble"))
for (m in 3:8){
  w<-ecalp(DATA2[,m], DATA2[,2])
  w<-as.data.frame(w)
  w$bin<-c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95)
  caltest<-lm(w~bin, data=w, subset=w>-1)
  ctest<-cor(DATA2[,m], DATA2[,2])
  Model_calibration_table_train$COR[m-2]<-ctest
  Model_calibration_table_train$calibration[m-2]<-caltest$coefficients[2]
  Model_calibration_table_train$bias[m-2]<-caltest$coefficients[1]
}
Presence_absence_thresholds_train<-optimal.thresholds(DATA2)
tht<-as.numeric(Presence_absence_thresholds_train[6,2:7])
Model_evaluation_table_train<-presence.absence.accuracy(DATA2, threshold=tht)
Model_evaluation_table_train$TSS<-Model_evaluation_table_train$sensitivity+Model_evaluation_table_train$specificity-1





#########################  PREDICTION TO TEST DATA SET ########################################

TEST$BRT_pred<-predict.gbm(BRT, TEST, n.trees=BRT$gbm.call$best.trees, type="response")
TEST$RF_pred<-predict(RF, TEST)
TEST$GLM_pred<-predict(GLM, type="response", newdata=TEST)
TEST$GAM_pred<-as.numeric(predict(GAM, type="response", newdata=TEST))
test2<-TEST
test2$SPECIES<-NULL
test2$Year<-as.numeric(test2$Year)
test2$season<-as.numeric(test2$season)
test2$COUNT<-NULL
TEST$MAX_pred<- predict(MAX, test2, factors=c('season', 'Year'))
rm(test2)

##### ENSEMBLE PREDICTION requires AUC values to calculate the weights for the ensemble predictionweights<-(Model_evaluation_table[,7]-0.5)/(sum(Model_evaluation_table[,7]-0.5))
DATA<-data.frame(TEST$LAT, TEST$SPECIES, TEST$BRT_pred, TEST$RF_pred, TEST$GLM_pred, TEST$GAM_pred, TEST$MAX_pred)
Model_evaluation_table<-presence.absence.accuracy(DATA)
weights<-(Model_evaluation_table[,7]-0.5)/(sum(Model_evaluation_table[,7]-0.5))
TEST$Ensemble_pred<-((TEST$BRT_pred*weights[1])+(TEST$RF_pred*weights[2])+ (TEST$GLM_pred*weights[3]) + (TEST$GAM_pred*weights[4]) + (TEST$MAX_pred*weights[5]))





##############  MODEL VALIDATION AND CALIBRATION OF TEST DATA FIT #############################
### NOTE: the optimal threshold and all threshold-dependent evaluation criteria were removed from the manuscript


DATA<-data.frame(TEST$LAT, TEST$SPECIES, TEST$BRT_pred, TEST$RF_pred, TEST$GLM_pred, TEST$GAM_pred, TEST$MAX_pred, TEST$Ensemble_pred)
Presence_absence_thresholds<-optimal.thresholds(DATA)
th<-as.numeric(Presence_absence_thresholds[6,2:7])
Model_evaluation_table<-presence.absence.accuracy(DATA, threshold=th)
Model_evaluation_table$TSS<-Model_evaluation_table$sensitivity+Model_evaluation_table$specificity-1
Model_calibration_table<-data.frame("Model"=c("BRT", "RF", "GLM", "GAM", "MAX", "Ensemble"))
for (m in 3:8){
  w<-ecalp(DATA[,m], DATA[,2])
  w<-as.data.frame(w)
  w$bin<-c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95)
  caltest<-lm(w~bin, data=w, subset=w>-1)
  ctest<-cor(DATA[,m], DATA[,2])
  Model_calibration_table$COR[m-2]<-ctest
  Model_calibration_table$calibration[m-2]<-caltest$coefficients[2]
  Model_calibration_table$bias[m-2]<-caltest$coefficients[1]
}





##############  SAVE OUTPUT TABLES ###########################################################

write.table(Model_evaluation_table_train,"Model_evaluation_table_train.csv", sep=",", row.names=F, quote=F)
write.table(Model_evaluation_table,"Model_evaluation_table.csv", sep=",", row.names=F, quote=F)
write.table(Presence_absence_thresholds_train,"Presence_absence_thresholds_train.csv", sep=",", row.names=F, quote=F)
write.table(Presence_absence_thresholds,"Presence_absence_thresholds.csv", sep=",", row.names=F, quote=F)
write.table(Model_calibration_table,"Model_calibration_table.csv", sep=",", row.names=F, quote=F)
write.table(Model_calibration_table_train,"Model_calibration_table_train.csv", sep=",", row.names=F, quote=F)
write.table(Model_calibration_table,"Model_calibration_table.csv", sep=",", row.names=F, quote=F)










##############################################################################################
#########################  CONSTRUCTING THE DENSITY MODELS ###################################
###################### WITH SEASON AND YEAR AS FACTOR PREDICTORS   ###########################
##############################################################################################



#########################  SETTING TRAINING AND TEST DATA  ###################################

t1<-subset(subset(SPEA, LAT>39.7), COUNT>0)
t2<-subset(subset(SPEA, LAT<38), COUNT>0)
TRAINa<-merge(t1, t2, all=T)
TRAINa<-subset(TRAINa, COUNT<500)                      		  #### remove extreme outliers
remove(t1, t2)



#########################  BOOSTED REGRESSION TREES ##########################################

BRTa <- gbm.step(data=TRAINa, gbm.x = 3:20, gbm.y = 1, family = "gaussian", tree.complexity = 8, learning.rate = 0.001, bag.fraction = 0.64)


#########################  RANDOM FOREST #####################################################

RFa<-randomForest(COUNT~LONG+LAT+effort+season+Year+BATIMETRY+front_distance+front_dens+gradient+DIST+SST+SSTgrad+SSTchange+SSTanomaly+CHL+CHLgrad+CHLchange+SSH, data=TRAINa, importance=T, replace=T, mtry=10, ntree=1500, na.action=na.omit)


#########################  GLM with NEGATIVE BINOMIAL ########################################
### NOTE: this model was removed from the manuscript

global<-glm.nb(COUNT~LONG+LAT+season+Year+offset(effort)+BATIMETRY+front_distance+front_dens+DIST+I(DIST^2)+SST+SSTgrad+SSTchange+SSTanomaly+CHLgrad+CHLchange+CHL+SSH, data=TRAINa, link = log)
NEGBIN<-stepAIC(global, trace = 0)



#########################  GLM with POISSON ##################################################

poisson<-glm(COUNT~LONG+LAT+season+Year+offset(effort)+BATIMETRY+front_distance+front_dens+DIST+I(DIST^2)+SST+SSTgrad+SSTchange+SSTanomaly+CHLgrad+CHLchange+CHL+SSH, data=TRAINa,  family=poisson(link=log))
POIS<-stepAIC(poisson, trace = 0)



#########################  GAM ###############################################################

GAMa<-gam(BALEARIC~Year+season+s(LONG,k=5)+s(LAT,k=5)+s(BATIMETRY, k=5)+s(front_distance,k=5)+s(front_dens,k=5)+s(gradient, k=5)+s(effort, k=5)+s(DIST, k=5)+s(SST, k=5)+s(SSH, k=5)+s(SSTgrad, k=5)+s(SSTchange, k=5)+s(SSTanomaly, k=5)+s(CHL, k=5)+s(CHLgrad, k=5)+s(CHLchange, k=5), data=TRAINa, family=poisson,optimizer="perf")



##############################################################################################
#################  ZERO INFLATED AND ZERO ALTERED MODELS #####################################
##############################################################################################
### NOTE: these models attempt to automatically model a distribution (0/1) and abundance component, and therefore use the full data set


#########################  ZERO-INFLATED POISSON #############################################
### variables for the presence part of the model are those selected for the distribution GLM
library(pscl)
ZIP<-zeroinfl(BALEARIC~LONG+LAT+BATIMETRY+season+Year+front_distance+front_dens+DIST+I(DIST^2)+SST+SSTgrad+SSTchange+SSTanomaly+CHLgrad+CHLchange+CHL+SSH | LAT+SST+season+Year+BATIMETRY+I(DIST^2)+DIST+front_dens, data=TRAIN)


#########################  ZERO-ALTERED OR HURDLE MODEL ######################################
### NOTE: this model was not included in the the manuscript

HURDLE<-hurdle(BALEARIC~LONG+LAT+BATIMETRY+season+Year+front_distance+front_dens+SST+CHL+SSH, data=TRAIN, dist='poisson', zero.dist='binomial', link='logit')







#########################  PREDICTION TO TRAINING DATA SET ####################################
##### TRAINING data are exclusively presence data, hence no presence/absence component required in prediction

TRAINa$BRT_pred<-predict.gbm(BRTa, TRAINa, n.trees=BRTa$gbm.call$best.trees, type="response")
TRAINa$RF_pred<-predict(RFa, TRAINa)
TRAINa$ZIP_pred<- predict(ZIP, type="response", newdata=TRAINa)
TRAINa$POIS_pred<- predict(POIS, type="response", newdata=TRAINa)

detach("package:pscl",unload=TRUE)						### package 'pscl' interferes with the 'gam' function in package 'mgcv' and must be detached
detach("package:gam",unload=TRUE) 						### package 'gam' interferes with the 'gam' function in package 'mgcv' and must be detached
TRAINa$GAM_pred<- as.numeric(predict(GAMa, type="response", newdata=TRAINa))


### The following will not be used further:
# TRAINa$NEGBIN_pred<-predict(NEGBIN, type="response", newdata=TRAINa) # removed during revision to make space for zero inflated model
# TRAINa$HURDLE_pred<- predict(HURDLE, type="response", newdata=TRAINa)



##### ENSEMBLE PREDICTION requires COR values to calculate the weights for the ensemble prediction:
DATAa2<-data.frame(TRAINa$COUNT, TRAINa$BRT_pred, TRAINa$RF_pred, TRAINa$ZIP_pred, TRAINa$POIS_pred, TRAINa$GAM_pred)
Model_calibration_table_abundance_train<-data.frame("Model"=c("BRT", "RF", "ZIP", "POIS", "GAM"))
for (m in 2:6){
  ctest<-cor(DATAa2[,m], DATAa2[,1])
  Model_calibration_table_abundance_train$COR[m-1]<-ctest
}
weights_t<-(Model_calibration_table_abundance_train[,2])/(sum(Model_calibration_table_abundance_train[,2]))

TRAINa$Ensemble_pred<-(TRAINa$BRT_pred*weights_t[1])+(TRAINa$RF_pred*weights_t[2])+ (TRAINa$ZIP_pred*weights_t[3]) + (TRAINa$POIS_pred*weights_t[4]) + (TRAINa$GAM_pred*weights_t[5])





##############  MODEL VALIDATION AND CALIBRATION OF TRAINING DATA FIT ##########################

DATAa2<-data.frame(TRAINa$COUNT, TRAINa$BRT_pred, TRAINa$RF_pred, TRAINa$ZIP_pred, TRAINa$POIS_pred, TRAINa$GAM_pred, TRAINa$Ensemble_pred)
Model_calibration_table_abundance_train<-data.frame("Model"=c("BRT", "RF", "ZIP", "POIS", "GAM", "Ensemble"))
for (m in 2:7){
  pred<-DATAa2[,m]
  obs<-DATAa2[,1]
  caltest<-lmodel2(log(obs+0.000001)~log(pred), nperm=100)			### regression is performed on log-transformed data, requiring a small addition to avoid -inf when data=0 
  ctest<-cor(DATAa2[,m], DATAa2[,1])
  Model_calibration_table_abundance_train$COR[m-1]<-ctest
  Model_calibration_table_abundance_train$R_sq[m-1]<-caltest$rsquare
  Model_calibration_table_abundance_train$calibration[m-1]<-caltest$regression.results[2,3]
  Model_calibration_table_abundance_train$bias[m-1]<-caltest$regression.results$Intercept[2]
}
write.table(Model_calibration_table_abundance_train,"Model_calibration_table_abundance_train.csv", sep=",", row.names=F, quote=F)





#########################  PREDICTION TO TEST DATA SET ########################################
##### TEST data include absences, hence the prediction is a combination of presence and abundance prediction #####
##### Alternatively, abundance could only be predicted to above-threshold presence predictions, with optimal thresholds taken from 'Presence_absence_thresholds  ##############################
##### 2=Sens=Spec, 3=Max(Sens+Spec), 4=MaxKappa, 5=MaxPCC, 6=PredPrev=obsPrev, 9=min(ROCdist) #


TEST$BRT_preda<-predict.gbm(BRTa, TEST, n.trees=BRTa$gbm.call$best.trees, type="response")*TEST$BRT_pred
TEST$RF_preda<-predict(RFa, TEST)*TEST$RF_pred
library(pscl)									### reload package 'pscl' (had to be detached for gam prediction)
TEST$ZIP_preda<-predict(ZIPois, type="response", newdata=TEST)    		### does not need to be multiplied as zero-inflated model takes zero's into account
TEST$POIS_preda<- predict(POIS, type="response", newdata=TEST)*TEST$GLM_pred

detach("package:pscl",unload=TRUE)						### package 'pscl' interferes with the 'gam' function in package 'mgcv' and must be detached
detach("package:gam",unload=TRUE) 						### package 'gam' interferes with the 'gam' function in package 'mgcv' and must be detached
TEST$GAM_preda<- as.numeric(predict(GAMa, type="response", newdata=TEST))*(TEST$GAM_pred)



### ENSEMBLE PREDICTION requires COR values to calculate the weights for the ensemble prediction:
DATAa<-data.frame(TEST$COUNT, TEST$BRT_preda, TEST$RF_preda, TEST$ZIP_preda, TEST$POIS_preda, TEST$GAM_preda)
Model_calibration_table_abundance<-data.frame("Model"=c("BRT", "RF", "ZIP", "POIS", "GAM"))
for (m in 2:6){
  ctest<-cor(DATAa[,m], DATAa[,1])
  Model_calibration_table_abundance$COR[m-1]<-ctest
}
weights<-(Model_calibration_table_abundance[,2])/(sum(Model_calibration_table_abundance[,2]))
TEST$Ensemble_preda<-(TEST$BRT_preda*weights[1])+(TEST$RF_preda*weights[2])+ (TEST$NEGBIN_preda*weights[3]) + (TEST$POIS_preda*weights[4]) + (TEST$GAM_preda*weights[5])






##############  MODEL VALIDATION AND CALIBRATION OF TEST DATA FIT #########################

DATAa<-data.frame(TEST$COUNT, TEST$BRT_preda, TEST$RF_preda, TEST$NEGBIN_preda, TEST$POIS_preda, TEST$GAM_preda, TEST$Ensemble_preda)
Model_calibration_table_abundance<-data.frame("Model"=c("BRT", "RF", "NEGBIN", "POIS", "GAM", "Ensemble"))
for (m in 2:7){
  pred<-DATAa[,m]
  obs<-DATAa[,1]
  caltest<-lmodel2(log(obs+0.000001)~log(pred), nperm=100)			### regression is performed on log-transformed data, requiring a small addition to avoid -inf when data=0 
  ctest<-cor(DATAa[,m], DATAa[,1])
  Model_calibration_table_abundance$COR[m-1]<-ctest
  Model_calibration_table_abundance$R_sq[m-1]<-caltest$rsquare
  Model_calibration_table_abundance$calibration[m-1]<-caltest$regression.results[2,3]
  Model_calibration_table_abundance$bias[m-1]<-caltest$regression.results$Intercept[2]
}
write.table(Model_calibration_table_abundance,"Model_calibration_table_abundance.csv", sep=",", row.names=F, quote=F)










###############################################################################################
##############  PREDICT DISTRIBUTION AND DENSITY TO ENTIRE STUDY AREA #########################
###############################################################################################



##############  LOAD ENVIRONMENTAL DATA FOR ALL SEASONS AND YEARS #############################

channel1<-odbcConnectAccess('SPEA_seabird_data.mdb')
FINAL_PREDICTION_EVERYWHERE <- sqlQuery(channel1, "SELECT * FROM PORTUGAL_all_years")
odbcClose(channel1)




##############  USE ABOVE MODELS TO PREDICT DISTRIBUTION AND DENSITY ##########################
############## Variable names and types must be identical to training data (except Maxent) ####

FINAL_PREDICTION_EVERYWHERE$MAX_pred<- predict(MAX, me, factors=c('season', 'Year'))	### Maxent prediction requires all numeric variables, hence is done first
FINAL_PREDICTION_EVERYWHERE$season<-as.factor(FINAL_PREDICTION_EVERYWHERE$season)	### After Maxent prediction, convert to factor
FINAL_PREDICTION_EVERYWHERE$Year<-as.factor(FINAL_PREDICTION_EVERYWHERE$Year)		### convert to factor
FINAL_PREDICTION_EVERYWHERE$BRT_pred<-predict.gbm(BRT, FINAL_PREDICTION_EVERYWHERE, n.trees=BRT$gbm.call$best.trees, type="response")
FINAL_PREDICTION_EVERYWHERE$RF_pred<-predict(RF, FINAL_PREDICTION_EVERYWHERE)
FINAL_PREDICTION_EVERYWHERE$GLM_pred<-predict(GLM, newdata=FINAL_PREDICTION_EVERYWHERE, type="response")
FINAL_PREDICTION_EVERYWHERE$GAM_pred<-as.numeric(predict(GAM, newdata=FINAL_PREDICTION_EVERYWHERE, type="response"))

##### ENSEMBLE PREDICTION BASED ON AUC-weighted averaging ###
f<-Model_evaluation_table[1:5,7]
weights<-(f-0.5)/(sum(f-0.5))
FINAL_PREDICTION_EVERYWHERE$Ensemble_pred<-((FINAL_PREDICTION_EVERYWHERE$BRT_pred*weights[1])+(FINAL_PREDICTION_EVERYWHERE$RF_pred*weights[2])+ (FINAL_PREDICTION_EVERYWHERE$GLM_pred*weights[3]) + (FINAL_PREDICTION_EVERYWHERE$GAM_pred*weights[4]) + (FINAL_PREDICTION_EVERYWHERE$MAX_pred*weights[5]))




##############  DENSITY PREDICTION by multiplying predicted abundance with predicted probability of presence ###############

FINAL_PREDICTION_EVERYWHERE$BRT_pred_ab<-FINAL_PREDICTION_EVERYWHERE$BRT_pred*predict.gbm(BRTa, FINAL_PREDICTION_EVERYWHERE, n.trees=BRT$gbm.call$best.trees, type="response")
FINAL_PREDICTION_EVERYWHERE$RF_pred_ab<-FINAL_PREDICTION_EVERYWHERE$RF_pred*predict(RFa, FINAL_PREDICTION_EVERYWHERE)
FINAL_PREDICTION_EVERYWHERE$ZIP_pred_ab<-predict(ZIP, type="response", newdata=FINAL_PREDICTION_EVERYWHERE)
FINAL_PREDICTION_EVERYWHERE$POIS_pred_ab<-FINAL_PREDICTION_EVERYWHERE$GLM_pred*predict(POIS, type="response", newdata=FINAL_PREDICTION_EVERYWHERE)
FINAL_PREDICTION_EVERYWHERE$GAM_pred_ab<-FINAL_PREDICTION_EVERYWHERE$GAM_pred*as.numeric(predict(GAMa, type="response", newdata=FINAL_PREDICTION_EVERYWHERE))



##### ENSEMBLE PREDICTION BASED ON COR-weighted averaging #######

weightsa<-(Model_calibration_table_abundance[,2])/sum(Model_calibration_table_abundance[,2])
FINAL_PREDICTION_EVERYWHERE$Ensemble_pred_ab<-((FINAL_PREDICTION_EVERYWHERE$BRT_pred_ab*weightsa[1])+(FINAL_PREDICTION_EVERYWHERE$RF_pred_ab*weightsa[2])+ (FINAL_PREDICTION_EVERYWHERE$ZIP_pred_ab*weightsa[3]) + (FINAL_PREDICTION_EVERYWHERE$POIS_pred_ab*weightsa[4]) + (FINAL_PREDICTION_EVERYWHERE$GAM_pred_ab*weightsa[5]))




##############  SAVE OUTPUT TABLES FOR USE IN GIS ############################################

write.dbf(FINAL_PREDICTION_EVERYWHERE, "EVERYWHERE_SPECIES_all_predictions.dbf", factor2char = TRUE, max_nchar = 254)







###############################################################################################
##############  USE ZONATION ALGORITHM TO IDENTIFY PRIORITY AREAS FOR CONSERVATION ############
###############################################################################################



##############  PREPARATION OF DATA FOR ZONATION SOFTWARE #####################################


library(SDMTools)
library(rgdal)
library(maptools)
library(colorspace)

setwd("C://Zonation")

for(m in 8:13){                                                            ### select the columns with predicted probabilities (here loop over 6 methods)
  spp_file<-data.frame()
  for (y in 2005:2009){
    yearf<-subset(FINAL_PREDICTION_EVERYWHERE, Year==y)
    
    Spring<-aggregate(yearf[,m]~LAT+LONG, data=yearf, subset=season==1, FUN=mean)
    Summer<-aggregate(yearf[,m]~LAT+LONG, data=yearf, subset=season==2, FUN=mean)
    Fall<-aggregate(yearf[,m]~LAT+LONG, data=yearf, subset=season==3, FUN=mean)
    
    ### CREATE ASCII FILES FROM DATA FRAMES
    
    dataframe2asc(Spring, filenames=sprintf("BASH_spring%i%s.asc",y,names(yearf)[m]), outdir='C://Zonation',gz = FALSE)
    dataframe2asc(Summer, filenames=sprintf("BASH_summer%i%s.asc",y,names(yearf)[m]), outdir='C://Zonation',gz = FALSE)
    dataframe2asc(Fall, filenames=sprintf("BASH_fall%i%s.asc",y,names(yearf)[m]), outdir='C://Zonation',gz = FALSE)
    
    ### CREATE SPECIES LIST FILE
    
    Season_weights<-c(1.0,1.0,1.0)  					### equal weights given to all three seasons
    alpha<-c(1.0,1.0,1.0)           					### no distribution smoothing used
    BQP<-c(1,1,1)
    buffer<-c(5,5,5)   							### buffer of surrounding 24 cells defined as affecting quality of target cell
    ABF<-c(1.0,1.0,1.0)           						### no additive benefit function used
    spec<-c(sprintf("BASH_spring%i%s.asc",y,names(yearf)[m]), sprintf("BASH_summer%i%s.asc",y,names(yearf)[m]), sprintf("BASH_fall%i%s.asc",y,names(yearf)[m]))
    spp<-data.frame(Season_weights, alpha, BQP, buffer, ABF, spec)
    spp_file<-rbind(spp_file, spp)
  }
  write.table(spp_file, sprintf("BASH_multiyear%s.spp",names(yearf)[m]), col.names=F, row.names=F, sep=" ", quote=F)
}



##############################################################################################################
##############################################################################################################
##############################################################################################################
# Run 'Zonation' for all six methods
# NOTE THAT THE "MODEL_run_settings.dat" is changed manually using Notepad

shell("C:/Zonation/BASH.bat", wait=T)

##############################################################################################################
##############################################################################################################
##############################################################################################################




##############  IMPORT Results Ascii file AND PLOT DATA FOR EACH METHOD #####################

setwd("C://Zonation")

labels<-c("BRT","RF","GLM","GAM","MaxEnt","Ensemble")
par(mfrow=c(3,2))
for(m in 8:13){
  yearf<-subset(FINAL_PREDICTION_EVERYWHERE, Year==2005)
  output<-aggregate(yearf[,m]~LAT+LONG, data=yearf, subset=season==1, FUN=mean)
  
  ### extract the data ###
  ### care needs to be taken to specify the files correctly! Zonation changes the file names depending on BLP settings
  ### here we use the files for the most recent run with BLP=0.1 and BQP=0
  
  
  ZONES<-readAsciiGrid(sprintf("%s.output.CAZ_EBLP50.rank.asc",names(yearf)[m]))  		### this reads it into a 'SpatialGridDataFrame' which is worthless for post-processing
  output$ZONES<-extract.data(output[,2:1],ZONES) 							### read out the data and write it into the original data frame
  
  ### draw a plot
  NS<-c(min(output$LAT),max(output$LAT))
  EW<-c(min(output$LONG),max(output$LONG))
  output$color<-level.colors(output$ZONES, at = do.breaks(c(0.75,1), 5),col.regions = colorRampPalette(sequential_hcl(5, h = 0, c = 0, l = c(100, 10), power = 1.5, fixup=T)))
  names(output)[3]<-labels[m-7]
  par(mar=c(1,1,0.1,0.1), plt= c(0.2, 0.8, 0.1, 0.9))
  plot(output$LAT~output$LONG, pch=15, cex=0.4, xlim=EW, ylim=NS, col=output$color, axes=T, ylab="", xlab="")
  text(-7,42.3,names(output)[3], cex=1.5)
  
}





##############################################################################################
#########################################  END ###############################################
###################### consider saving Workspace before closing R   ##########################
##############################################################################################

