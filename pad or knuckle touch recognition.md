

```R
library(caret)
library(ElemStatLearn)
library(rpart)
library(kernlab)
library(Hmisc)
library(e1071)
library(MASS)
library(ISLR)
library(boot)
library(quantmod) 
require(dplyr)
library(ggplot2)
library(psych)
library(polycor)
library(ISLR)
library(pROC)
library(ROCR)
library(Deducer)
library(randomForest)
library(tuneR)

```

    Loading required package: lattice
    Loading required package: ggplot2
    Loading required package: grid
    Loading required package: survival
    
    Attaching package: 'survival'
    
    The following object is masked from 'package:caret':
    
        cluster
    
    Loading required package: Formula
    
    Attaching package: 'Hmisc'
    
    The following objects are masked from 'package:base':
    
        format.pval, round.POSIXt, trunc.POSIXt, units
    
    
    Attaching package: 'e1071'
    
    The following object is masked from 'package:Hmisc':
    
        impute
    
    
    Attaching package: 'boot'
    
    The following object is masked from 'package:survival':
    
        aml
    
    The following object is masked from 'package:lattice':
    
        melanoma
    
    Loading required package: xts
    Loading required package: zoo
    
    Attaching package: 'zoo'
    
    The following objects are masked from 'package:base':
    
        as.Date, as.Date.numeric
    
    Loading required package: TTR
    Version 0.4-0 included new data defaults. See ?getSymbols.
    
    Attaching package: 'quantmod'
    
    The following object is masked from 'package:Hmisc':
    
        Lag
    
    Loading required package: dplyr
    
    Attaching package: 'dplyr'
    
    The following objects are masked from 'package:xts':
    
        first, last
    
    The following object is masked from 'package:MASS':
    
        select
    
    The following objects are masked from 'package:Hmisc':
    
        combine, src, summarize
    
    The following objects are masked from 'package:stats':
    
        filter, lag
    
    The following objects are masked from 'package:base':
    
        intersect, setdiff, setequal, union
    
    
    Attaching package: 'psych'
    
    The following object is masked from 'package:boot':
    
        logit
    
    The following object is masked from 'package:Hmisc':
    
        describe
    
    The following object is masked from 'package:kernlab':
    
        alpha
    
    The following object is masked from 'package:ggplot2':
    
        %+%
    
    Loading required package: mvtnorm
    Loading required package: sfsmisc
    
    Attaching package: 'sfsmisc'
    
    The following object is masked from 'package:dplyr':
    
        last
    
    The following object is masked from 'package:xts':
    
        last
    
    The following object is masked from 'package:Hmisc':
    
        errbar
    
    
    Attaching package: 'polycor'
    
    The following object is masked from 'package:psych':
    
        polyserial
    
    Type 'citation("pROC")' for a citation.
    
    Attaching package: 'pROC'
    
    The following objects are masked from 'package:stats':
    
        cov, smooth, var
    
    Loading required package: gplots
    
    Attaching package: 'gplots'
    
    The following object is masked from 'package:stats':
    
        lowess
    
    Loading required package: JGR
    Loading required package: rJava
    


    Error: package 'rJava' could not be loaded
    


    randomForest 4.6-12
    Type rfNews() to see new features/changes/bug fixes.
    
    Attaching package: 'randomForest'
    
    The following object is masked from 'package:psych':
    
        outlier
    
    The following object is masked from 'package:dplyr':
    
        combine
    
    The following object is masked from 'package:Hmisc':
    
        combine
    
    tuneR >= 1.0 has changed its Wave class definition.
    Use updateWave(object) to convert Wave objects saved with previous versions of tuneR.
    

## A BASIC FINGERSENSE CLASSIFIER:
### The task is to write a program that takes as input files describing containing information about a touch point and determine whether the touch was from a pad or knuckle touch. The program will classify approximately 10,000 touches. 


## DATASETS:
### The task_data zip file contains two additional zip files, train.zip and test.zip. Both train.zip and test.zip have the same directory structure.

## DIRECTORY STRUCTURE:

### The structure of train.zip and test.zip is as follows: 
### root
    ###  user_folder: [hand,table]-timestamp
        ### instance_folder: timestamp-[pad,knuckle]
            ### audio.wav
            ### touch.csv


### An instance represents data from a single finger tap. An instance contains information about the touch (x,y,touch major and minor axes, pressure, and orientation). Each tap may be from a pad or knuckle. Each instance folder represents a single instance; its label (pad/knuckle) is specified in directory name. Thetime stamp for each instance is guaranteed to be unique.Each user folder contains a collection of instance subfolders. A user folder represents data collected from a single user. The hand, table prefix on the user folder specifies whether the data for this user was collected when the user was holding a device in his/her hand, or resting the device on the table. The set of users in train and test are disjoint—no user in train is also in test. The instance folders in the test dataset do not have labels. The task is to generate labels for these test instances. The train.zip file should contain 20,659 training instances (10,255 knuckle, 10404 pad). The test.zip file should contain 10,528 test instances (5,263 knuckle, 5,265 pad).




```R
setwd("E:/Work/Rstudio/WorkDirectory/New")

```


```R
getthewavevalues<-function(sndObj){#data_as_list_train
        s1 <- sndObj@left
        s1 <- s1 / 2^(sndObj@bit -1)
        timeArray <- (0:(256-1)) / sndObj@samp.rate
        timeArray <- timeArray * 1000
        #plot(timeArray, s1, type='l', col='black', xlab='Time (ms)', ylab='Amplitude') 
        n <- length(s1)
        p <- fft(s1)
        nUniquePts <- ceiling((n+1)/2)
        p <- p[1:nUniquePts] #select just the first half since the second half 
        # is a mirror image of the first
        p <- abs(p)  #take the absolute value, or the magnitude 
        p <- p / n #scale by the number of points so that
        # the magnitude does not depend on the length 
        # of the signal or on its sampling frequency  
        p <- p^2  # square it to get the power 
        
        # multiply by two (see technical document for details)
        # odd nfft excludes Nyquist point
        if (n %% 2 > 0){
                p[2:length(p)] <- p[2:length(p)]*2 # we've got odd number of points fft
        } else {
                p[2: (length(p) -1)] <- p[2: (length(p) -1)]*2 # we've got even number of points fft
        }
        
        freqArray <- (0:(nUniquePts-1)) * (sndObj@samp.rate / n) #  create the frequency array 
        
        #plot(freqArray/1000, 10*log10(p), type='l', col='black', xlab='Frequency (kHz)', ylab='Power (dB)')
        rms_val <- sqrt(mean(s1^2))
        return(sqrt(sum(p)))
}
#getting the sound 
get.sound.test<-function(path="./test"){
        #Creating the list of *.wav files from train folder :
        list_of_files_test <- dir(path=path,pattern='.*[.]wav', recursive = T,full.names=TRUE)
        #Reading *.wav files into one list:
        data_as_list_test = lapply(list_of_files_test, function(x){readWave(file=x)})
        #making usable dataframe for train set:
        data_l=c()
        for (i in 1:length(data_as_list_test)) {
                data_l <- c(data_l,(getthewavevalues(data_as_list_test[[i]])))
        }
        return (data_l)
}


get.sound.train<-function(path="./train"){
        #Creating the list of train files:
        list_of_files_train <- dir(path=path,pattern='.*[.]wav', recursive = T,full.names=TRUE)
        #Reading train files into one list:
        data_as_list_train = lapply(list_of_files_train, function(x){readWave(file=x)})
        #making usable dataframe for train set:
        data_l=c()
        for (i in 1:length(data_as_list_train)) {
                data_l <- c(data_l,(getthewavevalues(data_as_list_train[[i]])))
        }
        return (data_l)
}

#MAKING TRAIN DATA FRAME:
get_train<-function(path){
        #Creating the list of train files:
        list_of_files_train <- dir(path=path,pattern='.*[.]csv', recursive = T,full.names=TRUE)
        
        #Reading train files into one list:
        data_as_list_train = lapply(list_of_files_train, function(x){read.csv(file=x,stringsAsFactors=F,na.strings=c(" ","None","Not Available"),header=T,sep=",")})
        #making usable dataframe for train set:
        df_train<-do.call(rbind, lapply(data_as_list_train, data.frame, stringsAsFactors=FALSE))
        
        #adding Furiee transformation variable to the train set
        df_train[,"fft"]<-get.sound.train("./train")
        
        #Adding user inforamtion and hand/table position:
        train_Positiocn_col<-gsub("./train/(.*)-(.*)-(.*)/touch.csv","\\1",list_of_files_train)
        df_train[,"Position"]<-train_Positiocn_col
        
        #Creating outcome column for train data frames:
        train_outcome_col<-gsub("./train.*-(.*).*/touch.csv","\\1",list_of_files_train)
        df_train[,"Outcome"]<-train_outcome_col
        df_train$Outcome<-as.factor(df_train$Outcome)
        
        return (df_train)
        
        
}

get_test<-function(path){
        #Creating the list of test files:
        list_of_files_test <- dir(path=path,pattern='.*[.]csv', recursive = T,full.names=TRUE)
        
        #Reading test files into one list
        data_as_list_test = lapply(list_of_files_test , function(x){read.csv(file=x,stringsAsFactors=F,na.strings=c(" ","None","Not Available"),header=T,sep=",")})
        
        #Making usable dataframe for train/test:
        df_test<-do.call(rbind, lapply( data_as_list_test, data.frame, stringsAsFactors=FALSE))
        
        #adding Furiee transformation variable to the test set
        df_test[,"fft"]<-get.sound.test("./test")
        
        #Adding user inforamtion and hand/table position:
        test_Position_col<-gsub("./test.*-(.*).*/touch.csv","\\1",list_of_files_test)
        df_test[,"Position"]<-test_Position_col
        return(df_test)
        
}
```


```R
train<-get_train("./train")
head(train,5)
```




<table>
<thead><tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>major</th><th scope=col>minor</th><th scope=col>pressure</th><th scope=col>orientation</th><th scope=col>fft</th><th scope=col>Position</th><th scope=col>Outcome</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>176</td><td>1647</td><td>6</td><td>6</td><td>0</td><td>-1</td><td>0.03354754</td><td>hand</td><td>pad</td></tr>
	<tr><th scope=row>2</th><td>152</td><td>1666</td><td>5</td><td>5</td><td>0</td><td>-1</td><td>0.03430865</td><td>hand</td><td>pad</td></tr>
	<tr><th scope=row>3</th><td>410</td><td>1703</td><td>5</td><td>5</td><td>0</td><td>-1</td><td>0.03601321</td><td>hand</td><td>pad</td></tr>
	<tr><th scope=row>4</th><td>439</td><td>1714</td><td>5</td><td>5</td><td>0</td><td>-1</td><td>0.04105755</td><td>hand</td><td>pad</td></tr>
	<tr><th scope=row>5</th><td>651</td><td>1692</td><td>6</td><td>6</td><td>0</td><td>-1</td><td>0.04558804</td><td>hand</td><td>pad</td></tr>
</tbody>
</table>





```R
head(train,4)
```




<table>
<thead><tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>major</th><th scope=col>minor</th><th scope=col>pressure</th><th scope=col>orientation</th><th scope=col>fft</th><th scope=col>Position</th><th scope=col>Outcome</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>176</td><td>1647</td><td>6</td><td>6</td><td>0</td><td>-1</td><td>0.03354754</td><td>hand</td><td>pad</td></tr>
	<tr><th scope=row>2</th><td>152</td><td>1666</td><td>5</td><td>5</td><td>0</td><td>-1</td><td>0.03430865</td><td>hand</td><td>pad</td></tr>
	<tr><th scope=row>3</th><td>410</td><td>1703</td><td>5</td><td>5</td><td>0</td><td>-1</td><td>0.03601321</td><td>hand</td><td>pad</td></tr>
	<tr><th scope=row>4</th><td>439</td><td>1714</td><td>5</td><td>5</td><td>0</td><td>-1</td><td>0.04105755</td><td>hand</td><td>pad</td></tr>
</tbody>
</table>





```R
test<-get_test("./test")


```


```R

head(test,5)
```




<table>
<thead><tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>major</th><th scope=col>minor</th><th scope=col>pressure</th><th scope=col>orientation</th><th scope=col>fft</th><th scope=col>Position</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>655</td><td>434</td><td>6</td><td>6</td><td>0</td><td>-1</td><td>0.04942666</td><td>20140213_135947/20140213_134857792</td></tr>
	<tr><th scope=row>2</th><td>176</td><td>411</td><td>6</td><td>6</td><td>0</td><td>-1</td><td>0.06784357</td><td>20140213_135947/20140213_134858346</td></tr>
	<tr><th scope=row>3</th><td>863</td><td>397</td><td>6</td><td>6</td><td>0</td><td>-1</td><td>0.05785467</td><td>20140213_135947/20140213_134859341</td></tr>
	<tr><th scope=row>4</th><td>686</td><td>684</td><td>4</td><td>3</td><td>0</td><td>-1</td><td>0.03590614</td><td>20140213_135947/20140213_134859819</td></tr>
	<tr><th scope=row>5</th><td>886</td><td>1038</td><td>3</td><td>3</td><td>0</td><td>-1</td><td>0.03594732</td><td>20140213_135947/20140213_134900662</td></tr>
</tbody>
</table>




### The train set will contain Outcome and Position columns. The Position colums added for consistency, to be sure that outcomes match corresponding instances and users .  Position values with hand/table values will be excluded from analisys as suggested in QeexoInterview.pdf


```R



```

### Exploratory data analysis:



### 1) No NAs are in data, so no imputing data needed:


```R
dim(train)
temp<-train[complete.cases(train),]
dim(temp)

```




<ol class=list-inline>
	<li>20659</li>
	<li>9</li>
</ol>







<ol class=list-inline>
	<li>20659</li>
	<li>9</li>
</ol>




### 2) Non-Zero variance variables: There two near-zero variance predictors: pressure and orientation 


```R
nzv <- nearZeroVar(train,saveMetrics= TRUE) 
nzv

```




<table>
<thead><tr><th></th><th scope=col>freqRatio</th><th scope=col>percentUnique</th><th scope=col>zeroVar</th><th scope=col>nzv</th></tr></thead>
<tbody>
	<tr><th scope=row>x</th><td>1.152174</td><td>5.218065</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>y</th><td>1</td><td>9.109831</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>major</th><td>1.861018</td><td>0.04840505</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>minor</th><td>1.949174</td><td>0.04356455</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>pressure</th><td>0</td><td>0.004840505</td><td>TRUE</td><td>TRUE</td></tr>
	<tr><th scope=row>orientation</th><td>0</td><td>0.004840505</td><td>TRUE</td><td>TRUE</td></tr>
	<tr><th scope=row>fft</th><td>1</td><td>99.94675</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>Position</th><td>1.376783</td><td>0.009681011</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>Outcome</th><td>1.014529</td><td>0.009681011</td><td>FALSE</td><td>FALSE</td></tr>
</tbody>
</table>




### Mean values and standard deviations for these valuaes 



```R
mean(train$pressure)
sd(train$pressure)
mean(train$orientation)
sd(train$orientation)

```




0






0






-1






0



### Standart deviations of the variables are both zeros and values of the variables do not change, so these variables will be excluded from predictions.



```R
train_set<-train[,c("x","y","major","minor","fft","Outcome")]
test_set<-test[,c("x","y","major","minor","fft")]

```

### 3) Two variables minor and major reveal correlation:


```R


```

### Correlation coefficient 0.85 confirms correlation between variables. The correlation will be  taken care of later in the models by introducint Principle Component Analisys into the models:



```R
cor(train_set[,-6])


```




<table>
<thead><tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>major</th><th scope=col>minor</th><th scope=col>fft</th></tr></thead>
<tbody>
	<tr><th scope=row>x</th><td> 1.00000000</td><td>-0.02727122</td><td>-0.02649693</td><td>-0.03879998</td><td>-0.12424263</td></tr>
	<tr><th scope=row>y</th><td>-0.02727122</td><td> 1.00000000</td><td>-0.05803024</td><td>-0.06204258</td><td>-0.13458393</td></tr>
	<tr><th scope=row>major</th><td>-0.02649693</td><td>-0.05803024</td><td> 1.00000000</td><td> 0.85893598</td><td>-0.20344000</td></tr>
	<tr><th scope=row>minor</th><td>-0.03879998</td><td>-0.06204258</td><td> 0.85893598</td><td> 1.00000000</td><td>-0.22654458</td></tr>
	<tr><th scope=row>fft</th><td>-0.1242426</td><td>-0.1345839</td><td>-0.2034400</td><td>-0.2265446</td><td> 1.0000000</td></tr>
</tbody>
</table>




### Creating train and test set:



```R

TrainIn<-createDataPartition(y=train_set[,dim(train_set)[2]],p=0.6,list=FALSE)
trains<-train_set[TrainIn,]
tests<-train_set[-TrainIn,]

dim(trains)
dim(tests)

```




<ol class=list-inline>
	<li>12396</li>
	<li>6</li>
</ol>







<ol class=list-inline>
	<li>8263</li>
	<li>6</li>
</ol>




### Histograms of the minor and major varaibles  have sckewed distribution, so the variables will be standardize by centering and scaling for the training set



```R
head(trains)
```




<table>
<thead><tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>major</th><th scope=col>minor</th><th scope=col>fft</th><th scope=col>Outcome</th></tr></thead>
<tbody>
	<tr><th scope=row>3</th><td>410</td><td>1703</td><td>5</td><td>5</td><td>0.03601321</td><td>pad</td></tr>
	<tr><th scope=row>4</th><td>439</td><td>1714</td><td>5</td><td>5</td><td>0.04105755</td><td>pad</td></tr>
	<tr><th scope=row>5</th><td>651</td><td>1692</td><td>6</td><td>6</td><td>0.04558804</td><td>pad</td></tr>
	<tr><th scope=row>6</th><td>651</td><td>1683</td><td>5</td><td>5</td><td>0.04913728</td><td>pad</td></tr>
	<tr><th scope=row>8</th><td>56</td><td>1427</td><td>6</td><td>5</td><td>0.04599574</td><td>pad</td></tr>
	<tr><th scope=row>10</th><td>421</td><td>1415</td><td>6</td><td>6</td><td>0.03948286</td><td>pad</td></tr>
</tbody>
</table>





```R
head(tests)
```




<table>
<thead><tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>major</th><th scope=col>minor</th><th scope=col>fft</th><th scope=col>Outcome</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>176</td><td>1647</td><td>6</td><td>6</td><td>0.03354754</td><td>pad</td></tr>
	<tr><th scope=row>2</th><td>152</td><td>1666</td><td>5</td><td>5</td><td>0.03430865</td><td>pad</td></tr>
	<tr><th scope=row>7</th><td>161</td><td>1430</td><td>6</td><td>5</td><td>0.0455487</td><td>pad</td></tr>
	<tr><th scope=row>9</th><td>457</td><td>1416</td><td>6</td><td>6</td><td>0.04089821</td><td>pad</td></tr>
	<tr><th scope=row>12</th><td>953</td><td>1633</td><td>5</td><td>5</td><td>0.04539964</td><td>pad</td></tr>
	<tr><th scope=row>16</th><td>659</td><td>1431</td><td>5</td><td>5</td><td>0.03438709</td><td>pad</td></tr>
</tbody>
</table>





```R

#Choosing a model:

#Logistic regression
ctrl <- trainControl(method = "cv", repeats = 10,classProbs = TRUE,number=10)

glmm<- train(Outcome ~ .,data = trains,method = "glm")
glmm_prob <- predict(glmm, newdata = tests,type = "prob")
glmm_train <- predict(glmm, newdata = trains)
glmm_pred <- predict(glmm, newdata = tests)


glm_trianCM <- confusionMatrix(trains$Outcome,glmm_train)
glm_CM <- confusionMatrix(tests$Outcome,glmm_pred)


```


```R
glm_trianCM
glm_CM
```




    Confusion Matrix and Statistics
    
              Reference
    Prediction knuckle  pad
       knuckle    5167  986
       pad        1680 4563
                                              
                   Accuracy : 0.7849          
                     95% CI : (0.7776, 0.7921)
        No Information Rate : 0.5524          
        P-Value [Acc > NIR] : < 2.2e-16       
                                              
                      Kappa : 0.5702          
     Mcnemar's Test P-Value : < 2.2e-16       
                                              
                Sensitivity : 0.7546          
                Specificity : 0.8223          
             Pos Pred Value : 0.8398          
             Neg Pred Value : 0.7309          
                 Prevalence : 0.5524          
             Detection Rate : 0.4168          
       Detection Prevalence : 0.4964          
          Balanced Accuracy : 0.7885          
                                              
           'Positive' Class : knuckle         
                                              






    Confusion Matrix and Statistics
    
              Reference
    Prediction knuckle  pad
       knuckle    3491  611
       pad        1094 3067
                                              
                   Accuracy : 0.7937          
                     95% CI : (0.7848, 0.8023)
        No Information Rate : 0.5549          
        P-Value [Acc > NIR] : < 2.2e-16       
                                              
                      Kappa : 0.5876          
     Mcnemar's Test P-Value : < 2.2e-16       
                                              
                Sensitivity : 0.7614          
                Specificity : 0.8339          
             Pos Pred Value : 0.8510          
             Neg Pred Value : 0.7371          
                 Prevalence : 0.5549          
             Detection Rate : 0.4225          
       Detection Prevalence : 0.4964          
          Balanced Accuracy : 0.7976          
                                              
           'Positive' Class : knuckle         
                                              




```R

```


```R

#Choosing a model:

#Logistic regression
ctrl <- trainControl(method = "cv", repeats = 20,classProbs = TRUE,number=20)

glmm<- train(Outcome ~ .,data = trains,method = "glm",trControl = ctrl, preProc = c("center", "scale","pca"))
glmm_prob <- predict(glmm, newdata = tests,type = "prob")
glmm_pred <- predict(glmm, newdata = tests)
glm_CM <- confusionMatrix(tests$Outcome,glmm_pred)

#Decision tree
DT<- train(Outcome ~ .,data = trains,method = "rpart",trControl = ctrl,preProc = c("center", "scale","pca"))
DT_prob <- predict(DT, newdata = tests, type = "prob")
DT_pred <- predict(DT, newdata = tests)
DT_CM <- confusionMatrix(tests$Outcome,DT_pred)

# lda
lda_fit<-train(Outcome ~ .,data = trains,method="lda",trControl = ctrl,preProc = c("center", "scale","pca"))
lda_prob<-predict(lda_fit, tests,type = "prob")
lda_pred<-predict(lda_fit, tests)
LDA_CM <- confusionMatrix(tests$Outcome,lda_pred)

#Random Forest
rft <- randomForest(Outcome ~ ., trains, ntree=100, norm.votes=FALSE)
rf_prob <- predict(rft, newdata = tests, type = "prob")
rf_pred<-predict(rft, tests)
rf_CM <- confusionMatrix(tests$Outcome,rf_pred)

#SVM
svmfit<-svm(Outcome~.,data=trains,cost=20,probability = T)
svm_tsp<-predict(svmfit,newdata=tests, probability = T)
svm_tsp_atr<-(attr(svm_tsp, "probabilities"))

svmfit_pf<-svm(Outcome~.,data=trains,cost=20,probability = F)
svm_ts_pf<-predict(svmfit_pf,newdata=tests, probability = F)
svm_ts_CM <- confusionMatrix(tests$Outcome,svm_ts_pf)




```


```R




```


```R

```


```R
df<-data.frame(glm_CM[3],DT_CM[3],LDA_CM[3],rf_CM[3],svm_ts_CM[3])
names<-c("Logistic regression","Decision tree","LDA","Random Forest","SVM")
names(df)<-names
df_t<-t(df)
dfn<-as.data.frame(df_t)
#dfn_s<-dfn[order(dfn$Accuracy),]
dfn[,1:2]
```




<table>
<thead><tr><th></th><th scope=col>Accuracy</th><th scope=col>Kappa</th></tr></thead>
<tbody>
	<tr><th scope=row>Logistic regression</th><td>0.7848239</td><td>0.5699678</td></tr>
	<tr><th scope=row>Decision tree</th><td>0.8158054</td><td>0.6322073</td></tr>
	<tr><th scope=row>LDA</th><td>0.7987414</td><td>0.5979063</td></tr>
	<tr><th scope=row>Random Forest</th><td>0.875832</td><td>0.7516728</td></tr>
	<tr><th scope=row>SVM</th><td>0.8404938</td><td>0.6808885</td></tr>
</tbody>
</table>




### The Accuracies for four chosen models gave us the best value for Random Forest  and the worst for LDA. In general we might have concluded that there are no drastic differences among all model, but let's build ROC curves:




```R
rocCurve_glm <- roc(response = tests$Outcome,predictor = glmm_prob[, "pad"],levels = rev(levels(tests$Outcome)))
rocCurve_lda <- roc(response = tests$Outcome,predictor = lda_prob[, "pad"],levels = rev(levels(tests$Outcome)))
rocCurve_DT <- roc(response = tests$Outcome,predictor = DT_prob[, "pad"],levels = rev(levels(tests$Outcome)))
rocCurve_rft <- roc(response = tests$Outcome,predictor = rf_prob[, "pad"],levels = rev(levels(tests$Outcome)))
rocCurve_svm<-roc(response = tests$Outcome,predictor = svm_tsp_atr[, "pad"],levels = rev(levels(tests$Outcome)))


```


```R
plot(rocCurve_glm)
plot(rocCurve_lda, add=TRUE, col='red')
plot(rocCurve_DT, add=TRUE, col='green')
plot(rocCurve_rft, add=TRUE, col='pink')
plot(rocCurve_svm, add=TRUE, col='blue')

```




    
    Call:
    roc.default(response = tests$Outcome, predictor = glmm_prob[,     "pad"], levels = rev(levels(tests$Outcome)))
    
    Data: glmm_prob[, "pad"] in 4161 controls (tests$Outcome pad) > 4102 cases (tests$Outcome knuckle).
    Area under the curve: 0.8952






    
    Call:
    roc.default(response = tests$Outcome, predictor = lda_prob[,     "pad"], levels = rev(levels(tests$Outcome)))
    
    Data: lda_prob[, "pad"] in 4161 controls (tests$Outcome pad) > 4102 cases (tests$Outcome knuckle).
    Area under the curve: 0.8946






    
    Call:
    roc.default(response = tests$Outcome, predictor = DT_prob[, "pad"],     levels = rev(levels(tests$Outcome)))
    
    Data: DT_prob[, "pad"] in 4161 controls (tests$Outcome pad) > 4102 cases (tests$Outcome knuckle).
    Area under the curve: 0.8166






    
    Call:
    roc.default(response = tests$Outcome, predictor = rf_prob[, "pad"],     levels = rev(levels(tests$Outcome)))
    
    Data: rf_prob[, "pad"] in 4161 controls (tests$Outcome pad) > 4102 cases (tests$Outcome knuckle).
    Area under the curve: 0.952






    
    Call:
    roc.default(response = tests$Outcome, predictor = svm_tsp_atr[,     "pad"], levels = rev(levels(tests$Outcome)))
    
    Data: svm_tsp_atr[, "pad"] in 4161 controls (tests$Outcome pad) > 4102 cases (tests$Outcome knuckle).
    Area under the curve: 0.9314




![svg](output_40_5.svg%2Bxml)


### The models with the largest area under the ROC curves (for Random Forest and SVM) will be choosen for our test prediction. As we can see the Logistic regression gave us the lowest AUC value.



```R
#dfn<-c("Random Forest", "Decision tree","Logistic regression","LDA")
AUC<-c(rocCurve_glm$auc,rocCurve_DT$auc,rocCurve_lda$auc,rocCurve_rft$auc,rocCurve_svm$auc)
dfn[,"AUC"]<-AUC
#dfnames<-c("Model","AUC","Accuracy")
FinalData<-dfn[,c("Accuracy","AUC")]
FinalData_o<-FinalData[order(-FinalData$Accuracy),]
FinalData_o



```




<table>
<thead><tr><th></th><th scope=col>Accuracy</th><th scope=col>AUC</th></tr></thead>
<tbody>
	<tr><th scope=row>Random Forest</th><td>0.875832</td><td>0.9520439</td></tr>
	<tr><th scope=row>SVM</th><td>0.8404938</td><td>0.9314374</td></tr>
	<tr><th scope=row>Decision tree</th><td>0.8158054</td><td>0.8166325</td></tr>
	<tr><th scope=row>LDA</th><td>0.7987414</td><td>0.8946327</td></tr>
	<tr><th scope=row>Logistic regression</th><td>0.7848239</td><td>0.8952481</td></tr>
</tbody>
</table>




### Conclusions: Based on AUC and Accuracy values we can conclude that SVM and Random Forest are the best chose for predicting our Outcome.


```R

```


```R

```


```R

```


```R

```


```R


```


```R


```


```R

```


```R

```


```R

```
