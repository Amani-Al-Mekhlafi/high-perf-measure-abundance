##### generate random data with the same missing pattern usining either Monte-Carlo or permutation test######### 
generate.random.mis <- function(data,method="Monte Carlo"){
## dattable is the data frame of biomarker candidates in columns, cases in rows, the last column has the state of disease   ##
## method is the method  to generate random data set either "MonteCarlo" or "Permutation"                                             ##
    
  labs <- data[,ncol(data)]
  dat<- data[,-ncol(data)]
  n<- nrow(data)
  m<- ncol(data[,-ncol(data)])
  
  
  if (method== "Monte Carlo"){
    
    scores <- matrix(runif(n*m), nrow=n, ncol=m)
    
    ind<-  unique(which(is.na(dat), arr.ind=TRUE))
    
    if (length(ind)>0){
      for(i in 1: nrow(ind)){
        scores[ind[i,1], ind[i,2]]<- NA
      }
    }
    
    return(data.frame(scores,labs))
  }
  
  if (method== "Permutation"){
    
    rand.labs<- sample(labs,size= n,replace=F)
    return(data.frame(dat,rand.labs))
  } 
  else{
    stop("the method of simulation should be either 'Monte Carlo' or 'Permutation' ")
  }

}


###############################################################################
###### For Misclassification rate as performance measurement ##################
####  only for multi-class problems, minimum 3 classes and maximum 10 #########
##############################################################################
Outlier.index<-function(x){

  quantiles<-quantile(x,na.rm=T)
  lowerq<-quantiles[2]
  upperq<-quantiles[4]
  iqr<-IQR(x,na.rm=T)
  extreme.upper<-upperq+1.5*iqr
  extreme.lower<-lowerq-1.5*iqr
  
  outliers<-NULL
  outliers<-c(which(x>extreme.upper),which(x<extreme.lower))
  return(outliers)
}

compute.mean<-function(x,labs){

  outliers<-NULL
  level<-levels(as.factor(labs))
  means<-rep(0,length(level))
  for(i in 1:length(level)){
    objs<-which(labs==level[i])
    outliers.class<-Outlier.index(x[objs])
    
    if(length(outliers.class)>0){
      means[i]<-mean(x[objs[-outliers.class]],na.rm=T)
      outliers<-c(outliers,objs[outliers.class])
    }else{
      means[i]<-mean(x[objs],na.rm=T)
    }
    
  }
  
  return(list(mean=means,outliers=outliers))
}

order.multiclass<-function(scores, labs){

  labs<-as.factor(labs)
  sort.scores<-sort(scores, index.return=T)
  
  scores<-sort.scores$x
  index<-sort.scores$ix
  labs<-labs[index]
  lv.lab<-levels(labs)
  res<-compute.mean(scores,labs)
  
  mean_lv<- res$mean
  names(mean_lv)<- lv.lab
  mean.lv<- sort(mean_lv)
  index.lv<-order(res$mean)
  
  outlier.index<-NULL
  if(!is.null(res$outliers)){
    outliers<-res$outliers
    scores<-scores[-outliers]
    labs<-labs[-outliers]
  }
  
  return(list(score=scores, lab=labs,level=lv.lab[index.lv], means= mean.lv))
}

check.similarity<- function(scores,labels){

  df<- data.frame(scores,labels)
  lv<- levels(as.factor(labels))
  logic<- NULL
  
  for (i in 1: length(lv)){
    first.clas<- subset(df,df[,ncol(df)]== lv[i])
    rest.class<- subset(df,df[,ncol(df)]!= lv[i])
    
    a<- any(first.clas$scores %in% rest.class$scores)
    logic[i]<- a
  }
  return(any(logic))
}

class.indx<- function(new.labs,lv){
  if(length(lv)<=2 |length(lv)>=11 ){
    stop("the number of classes should be more than 2 and less than 11")}
  n.cutoff<- length(lv)-1

  class.index<- rep(list(NULL), n.cutoff)
  
  lab.index<- which(new.labs== lv[2])
  lab.index1<- lab.index+1
  
  class.index[[1]]<- setdiff(lab.index, c(1,lab.index1))

  for(i in 2:n.cutoff){
    lab.index<-which(new.labs==lv[i+1])
    
    c.min.index<-class.index[[i-1]][1]
    large.lab.index<-which(lab.index>c.min.index)
    
    lab.index<-lab.index[large.lab.index]
    lab.index1<-lab.index+1
    
    class.index[[i]]<-setdiff(lab.index,lab.index1)
    
  }
  return(class.index)
}

add.noise<- function(scores){
 
  ## deal with NA
  scores[which(is.na(scores)== TRUE)]<- 0
  diff.scores<- diff(scores)
  ## deal with 0
  if(length(which(diff.scores==0))> 0){
    diff.score<- diff.scores[-c(which(diff.scores==0))] 
  }else{
    diff.score<- diff.scores
  }
  
  min.diff<- min(abs(diff.score))/2
  
  noised.scores<- NULL
  for(i in 1:length(scores)){
    x<- scores[i]+ runif(1,-min.diff,min.diff)
    noised.scores[i]<- x
  }
  
  return(noised.scores)
}

class.index.list<- function(scores,labels){

  if(check.similarity(scores,labels)){
    
    last.scores<- add.noise(scores)
    new.order<- order.multiclass(last.scores, labels)
    Index<- class.indx(new.order$lab,new.order$level)
    
  }else{
    new.order<- order.multiclass(scores, labels)
    Index<- class.indx(new.order$lab,new.order$level)
  }
  
  return(Index)
}

possible.cutpoint<- function(index.mat){

  n.cutoff<- length(index.mat)
  if( n.cutoff <2 |  n.cutoff>9){stop("the number of classes should be more than 2 and less than 11")}
  if(n.cutoff== 2){return(as.matrix(do.call(data.table::CJ, index.mat)[V1 < V2 ]))}
  if(n.cutoff== 3){return(as.matrix(do.call(data.table::CJ, index.mat)[V1 < V2 & V1 < V3 & V2 < V3]))}
  if(n.cutoff== 4){return(as.matrix(do.call(data.table::CJ, index.mat)[V1<V2 & V1<V3 & V1<V4 & V2<V3 & V2<V4 &V3<V4]))}
  if(n.cutoff== 5){return(as.matrix(do.call(data.table::CJ, index.mat)[V1<V2 &V1<V3 &V1<V4 &V1<V5 &V2<V3 &V2<V4 &V2<V5 &V3<V4 &V3<V5 &V4<V5]))}
  if(n.cutoff== 6){return(as.matrix(do.call(data.table::CJ, index.mat)[V1<V2&V1<V3&V1<V4&V1<V5&V1<V6&V2<V3&V2<V4&V2<V5&V2<V6&V3<V4&V3<V5&V3<V6&V4<V5&V4<V6&V5<V6]))}
  if(n.cutoff== 7){return(as.matrix(do.call(data.table::CJ, index.mat)[V1<V2&V1<V3&V1<V4&V1<V5&V1<V6&V1<V7&V2<V3&V2<V4&V2<V5&V2<V6&V2<V7&V3<V4&V3<V5&V3<V6&V3<V7&V4<V5&V4<V6&V4<V7&V5<V6&V5<V7&V6<V7]))}
  if(n.cutoff== 8){return(as.matrix(do.call(data.table::CJ, index.mat)[V1<V2&V1<V3&V1<V4&V1<V5&V1<V6&V1<V7&V1<V8&V2<V3&V2<V4&V2<V5&V2<V6&V2<V7&V2<V8&V3<V4&V3<V5&V3<V6&V3<V7&V3<V8&V4<V5&V4<V6&V4<V7&V4<V8&V5<V6&V5<V7&V5<V8&V6<V7&V6<V8&V7<V8]))}
  if(n.cutoff== 9){return(as.matrix(do.call(data.table::CJ, index.mat)[V1<V2&V1<V3&V1<V4&V1<V5&V1<V6&V1<V7&V1<V8&V1<V9&V2<V3&V2<V4&V2<V5&V2<V6&V2<V7&V2<V8&V2<V9&V3<V4&V3<V5&V3<V6&V3<V7&V3<V8&V3<V9&V4<V5&V4<V6&V4<V7&V4<V8&V4<V9&V5<V6&V5<V7&V5<V8&V5<V9&V6<V7&V6<V8&V6<V9&V7<V8&V7<V9&V8<V9]))}
  }

calculate.props<- function(cutpoint.matrix, n.scores){

  n.cutoff<- ncol(cutpoint.matrix)
  props<- matrix(0, nrow= nrow(cutpoint.matrix), ncol= n.cutoff+1)
  props[,1]<- cutpoint.matrix[,1]- 1
  
  
  for(i in 2:n.cutoff){
    
    props[,i]<- cutpoint.matrix[,i]- cutpoint.matrix[,i-1]
    
    
  }
  
  props[,n.cutoff+1]<- (n.scores - cutpoint.matrix[,n.cutoff])+1
  
  return(props)
}

calculate.false<- function(scores, labs, props){

  order.scores<- order(scores)
  order.labs<- labs[order.scores]
  
  false.labs<- NULL
  for (i in 1:nrow(props)){
    
    split.labs<- BBmisc::chunk(order.labs, props = c(props[i,]))
    
    false.sum<- c(1:length(levels(labs)))
    for(j in 1: ncol(props)){
      tab<-table(split.labs[[j]])
      
      a<- max(tab)
      b<- which(tab==a)
      false.sum[j]<- sum(tab)- tab[b]
    }
    
    false.labs[i]<- sum(false.sum)
  }
  
  return (false.labs)
}

minm.false<-function(scores, labs, props){

  xx<- suppressWarnings({calculate.false(scores, labs, props)})
  min.xx<- min(xx)
  minm.xx<- which(xx== min.xx)

  return (list(n.false= min.xx, cut.score= minm.xx))
}

misClass.Rate<- function(scores, labels){

  class.index<- class.index.list(scores,labels)
  cutpoint.matrix<- possible.cutpoint (class.index)
  props<- calculate.props(cutpoint.matrix, n.scores= length(scores))
  
  min.false<- minm.false(scores,labels, props) 
  mis.rate<- min.false$n.false/length(scores)

  return (mis.rate)

}


#################################################################################
#####  calculate the performance measurements    ##################
#############################################################################
aac<- function(x,labs, pos.class=NULL){
  labs<- as.factor(labs)
  lvs<- levels(labs)
  if (is.null(pos.class)){
    pos.class <- lvs[1]
  }else{
    pos.class<- pos.class
  }
  
  dat<- data.frame(x,labs)
  aac<- cost.curve(data=dat, attrs.no=1, pos.class)
    c<- dev.off()
  
  return(aac)
}

AAC <- function(scores,labs, pos.class){
  return(aac(scores, labs, pos.class))
}

HUM.perf<- function(x,labs){
  data<- data.frame(x,labs)
  indexClass=ncol(data)
  label=levels(data[,indexClass])
  indexLabs= label[1:length(label)]
  
  a<- CalculateHUM_seq(data,indexF=1,indexClass=ncol(data),indexLabel= indexLabs)

  return(a)
}

HUM <- function(scores,labs){
  return(HUM.perf(scores, labs)$HUM)
}

entropy <- function(scores, labs){
  return(cutIndex(scores, labs)[2])
}

mauc<- function(x, labs){
  
  a<- pROC:: multiclass.roc(labs, x)$auc[1]
  
  if (a>= 0.5){
    auc<- a
  }else{
    auc<- 1-a
  }
  
  return(a)
}
mAUC <- function(scores,labs){
  return(multiclass.roc(labs,scores)$auc[1])

}

misClassRate <- function(scores,labs){
  return(misClass.Rate(scores, labs))
}

perf.Measure<- function(data, method ,pos.class=NULL){
  data[,ncol(data)]<- as.factor(data[,ncol(data)])
  switch(method, 
         mAUC={mauc<- apply(data[,-ncol(data)],2, mAUC, data[,ncol(data)])
         vals<- mauc
         
         }, 
         entropy<- {ent<- apply(data[,-ncol(data)],2, cutIndex, data[,ncol(data)])
         if(is.list(ent)){
           vals<- unlist(sapply(ent, "[[", 2))}
         if(is.matrix(ent)){
           vals<- ent[2,] }
         },
         AAC={ vals<- apply(data[,-ncol(data)], 2, aac, data[,ncol(data)],  pos.class)    
         }, 
         HUM={hum<-apply(data[,-ncol(data)], 2, HUM.perf, data[,ncol(data)])
         vals<- unlist(sapply(hum, "[[", 1))
         },
         misClassRate={misR<-apply(data[,-ncol(data)], 2, misClassRate, data[,ncol(data)])
         vals<- misR
         }
         
  )
  return(vals)
  
}

### This function generated no.simulation "biomarkers" with random values
### that have no association with classes.

performance.measure.sim <- function(pfm,n,p.classes,no.simulations=1000, pos.class=NULL){
### pfm is the performance measure to be used for the simulated biomarkers.
### n is the sample size.
### p.classes is the distribution of the classes (giben in proportions)
## pos.class is the positive class, it is needed when the Area above the cost curve is measured, however in case it is written NULL the first class, ordered alphabetically, will be chosen.                                                                                           ##
  
  ### Generate the class labels.
  props <- p.classes/sum(p.classes)
  diagn <- rep(1,n)
  ind <- round(n*props[1])
  for (i in 2:length(props)){
    ind.start <- ind + 1
    ind <- round(n*sum(props[1:i]))
    diagn[ind.start:ind] <- i
  }
  
  diagn <- as.factor(diagn)

  
  ### Generate no.simulations "biomarkers" with random values.
  
  x <- matrix(runif(n*no.simulations),nrow=n)
  
  ### Compute the performance measure values for the generated "biomarkers".
  if(pfm== "AAC"){
    pfms <- sapply(data.frame(x),pfm,labs=diagn,pos.class= pos.class)
  }else{
    pfms <- sapply(data.frame(x),pfm,labs=diagn)
  }
  return(pfms)
}



plotLines<-function(x, y, xlab=NULL, ylab=NULL, main=NULL, titles=NULL,  shape="vh", opacity=NULL) {
  
  if(is.vector(y)) y<-matrix(y, nrow=length(x))
  
  if(is.null(titles)) titles<-paste(colnames(y))
  if(is.null(opacity)) opacity<-rep(1, ncol(y))
  
  
  if(length(x)!=nrow(y) | ncol(y)!=length(titles) | ncol(y)!=length(opacity))
    stop("The length of 'x', the row number of 'y' must be identical, as well as
         the column number of 'y', the length of 'titles' and 'opacity' !")
  
  
  
  xa <- list(
    title = xlab,
    titlefont = t
  )
  ya <- list(
    title = ylab,
    titlefont = t
  )
  
  
  p<-plotly::plot_ly(x=x) %>% plotly::layout(title=main,xaxis = xa, yaxis = ya)
  
  
  for(i in 1:ncol(y)) {
    #p<-plotly::add_lines(p, y=y[, i], line=list(shape=shape), opacity=opacity[i], name=titles[i])
    p<-plotly::add_lines(p, y=y[, i], line=list(shape=shape, width =6), opacity=opacity[i], name=titles[i])
  }
  p
}

 
performance<- function(dattable,imput.method, pfM.method, 
  is.positive, no.simulations,pos.class, corrected.method){

## Computes the performance of the biomarkers candidates 
## with the p-value that calculates by simulations and corrected by         ##
## Benjamini and Hochberg method or by Holm-Bonferroni                                                                            ##
##                                                                                                                                ##
## dattable              data frame of biomarker candidates in columns, cases in rows, the last column has the state of disease   ##
## imput.method          how to impute the missing values either "median" or "random"                                             ##
## pfM.method            performance measurement either "entropy","mAUC","AAC","HUM","misClassRate"                               ##
## no.simulations        number of simulation to calculate p-values                                                               ##
## pos.class             positive class                                                                                           ##
## is.positive           relation between the value of the performance measurement and the quality of the marker.                 ##
## corrected.method      corrected the p-values either by controlling FWER or FDR                                                 ## 
####################################################################################################################################
 
dattable[,ncol(dattable)]<- as.factor(dattable[,ncol(dattable)])
  
  ## 1. Missing Imputation
  if (sum(is.na(dattable))> 0){
    ## Real dattable
    indx <- unique(which(is.na(dattable), arr.ind=TRUE)[,2])
    misMat<- matrix(0, nrow(dattable), length(indx))
    colnames(misMat)<- colnames(dattable)[indx]
    
    
    if(imput.method== "median"){
      for(i in 1:length(indx)){
        a<- Hmisc::impute(dattable[,indx[i]])
        
        misMat[,i]<- a
      }
      
      dattable[,colnames(dattable)[indx]]<- misMat
      
    }
    ###################################
    if (imput.method== "random"){
      for(i in 1:length(indx)){
        a<- Hmisc::impute(dattable[,indx[i]], "random")
        misMat[,i]<- a
      }
      
      dattable[,colnames(dattable)[indx]]<- misMat
    }
    
  }

  ## 2. performance measurment
  perform<- perf.Measure(dattable, method= pfM.method ,pos.class) 

  ## 3. simulation
  labs<- dattable[,ncol(dattable)]
  tab<-prop.table(table(labs))
  pfms<-performance.measure.sim (pfm= pfM.method,nrow(dattable),tab,no.simulations= no.simulations) 

  ## 4. p values & confidence interval 
  if(is.positive==TRUE){

    pfms_sum<-sapply(perform, function(x) sum(pfms>=x))
    
  }
  
  if(is.positive==FALSE){

    pfms_sum<-sapply(perform, function(x) sum(pfms<=x))
  }

    p_value<- NULL
  for (i in 1: length(pfms_sum)){
    probability<- pfms_sum[i]/no.simulations
    p_value[i]<- probability
  }
  
  ## 5.Corrected p values
  
  if(corrected.method== "FWER"){
    correctd_P.value<- p_value* (ncol(dattable)-1)
  }
  
  if(corrected.method== "FDR"){
    correctd_P.value<- p.adjust(p_value, method = "BH", n = length(p_value))
  }
  
  correctd_P.value<- round(correctd_P.value,4)
  correctd_P.value<- ifelse(correctd_P.value  < 0.001, "< 0.001", correctd_P.value )
  correctd_P.value<- ifelse(correctd_P.value > 0.05,  "1", correctd_P.value )
  
  p_value<- round(p_value, 4)
  p_value<- ifelse(p_value < 0.001, "< 0.001", p_value )
  p_value<- ifelse(p_value > 0.05,  "> 0.05", p_value)

  d<- data.frame(colnames(dattable[,-ncol(dattable)]),perform,p_value,correctd_P.value)
  colnames(d)<- c("Biomarker candidates", paste(pfM.method), "p value", paste("corrected p value",corrected.method))
  rownames(d)<-NULL
  
  return(d)
  
}

hipermab<- function(dattable,random.simulation, imput.method, pfM.method,no.simulations,pos.class,Con.Interval,is.positive=TRUE){
                                                                                                                          
## dattable              data frame of biomarker candidates in columns, cases in rows, the last column has the state of disease   ##
## random.simulation     how to generate random data set either "MonteCarlo" or "Permutation"                                     ##
## imput.method          how to impute the missing values either "median" or "random"                                             ##
## pfM.method            performance measurement either "entropy","mAUC","AAC","HUM","misClassRate"                               ##
## no.simulations        number of simulation to calculate p-values                                                               ##
## pos.class             positive class                                                                                           ##
## Con.Interval          Confidence Interval                                                                                      ##
## is.positive           relation between the value of the performance measurement and the quality of the marker.                 ##
#################################################################################################################################### 
 
   labs<- as.factor(dattable[,ncol(dattable)])
  lvs<- levels(labs)
  
  if (is.null(pos.class)){
    pos.class <- lvs[1]
  }
  
  ## 1. Generating Random dattable
  rand.dat<- generate.random.mis(dattable,random.simulation)
  
  
  ## 2. Imputation
  if (sum(is.na(dattable))> 0){
    ## Real dattable
    indx <- unique(which(is.na(dattable), arr.ind=TRUE)[,2])
    misMat<- matrix(0, nrow(dattable), length(indx))
    colnames(misMat)<- colnames(dattable)[indx]
    
    ## Random dattable
    indx.rand <- unique(which(is.na(rand.dat), arr.ind=TRUE)[,2])
    misMat.rand<- matrix(0, nrow(rand.dat), length(indx.rand))
    colnames(misMat.rand)<- colnames(rand.dat)[indx.rand]
    
    if(imput.method== "median"){
      for(i in 1:length(indx)){
        a<- Hmisc::impute(dattable[,indx[i]])
        
        misMat[,i]<- a
      }
      for(i in 1:length(indx.rand)){
        a<- Hmisc::impute(rand.dat[,indx.rand[i]])
        misMat.rand[,i]<- a
      }
      
      rand.dat[,colnames(rand.dat)[indx.rand]]<- misMat.rand
      dattable[,colnames(dattable)[indx]]<- misMat
      
    }
    ###################################
    if (imput.method== "random"){
      for(i in 1:length(indx)){
        a<- Hmisc::impute(dattable[,indx[i]], "random")
        misMat[,i]<- a
      }
      for(i in 1:length(indx.rand)){
        a<- Hmisc::impute(rand.dat[,indx.rand[i]], "random")
        misMat.rand[,i]<- a
      }
      
      rand.dat[,colnames(rand.dat)[indx.rand]]<- misMat.rand
      dattable[,colnames(dattable)[indx]]<- misMat
    }
    
  }
  
  
  perform<- perf.Measure(dattable, method= pfM.method ,pos.class)
  perform.rand<- perf.Measure(rand.dat, method= pfM.method ,pos.class)
  
  ## 4. simulation
  tab<-prop.table(table(labs))
  pfms<-performance.measure.sim (pfM.method,nrow(dattable),tab,no.simulations= no.simulations)
  
  ## 5.  seqs
  all.ent<- c(min(perform), max(perform),min(perform.rand), max(perform.rand),min(pfms), max(pfms))
  seqs<- seq(min(all.ent), max(all.ent), 0.01)
  
  if(is.positive==TRUE){
    RealData<-sapply(seqs, function(x) sum(perform>x))
    RandomData<-sapply(seqs, function(x) sum(perform.rand>x))
    pfms_sum<-sapply(seqs, function(x) sum(pfms>x))
    
  }
  
  if(is.positive==FALSE){
    RealData<-sapply(seqs, function(x) sum(perform<x))
    RandomData<-sapply(seqs, function(x) sum(perform.rand<x))
    pfms_sum<-sapply(seqs, function(x) sum(pfms<x))
  }
  
  ##6.Con.Interval
  conf.Interval<- NULL
  
  for (i in 1: length(pfms_sum)){
    probability<- pfms_sum[i]/no.simulations
    CI<- qbinom(Con.Interval, ncol(dattable[,-ncol(dattable)]),probability )
    conf.Interval[i]<- CI
    
  }
  
  df<- data.frame(RealData, RandomData, conf.Interval)
  colnames(df)<- c("Real Data", "Random Data", "confidence Interval")
  rownames(df)<- seqs
  return(list(Performance= data.frame(perform),
              ByNumbers= df,
              ByPlot= plotLines(seqs, df,xlab=paste(pfM.method),  ylab="Number of Features", main="HiPerMAb Curves", titles=NULL, shape="vh")))
  
}

####################################################################################################################################
## Computes the probability that the no. of biomarker candidates with a good performance measure exceeds the (1-alpha) CI,        ##
## assuming we have m biomarker candidates out of which k have a p-value of p.0 (or lower).                                       ##
##                                                                                                                                ##
## Arguments                                                                                                                      ##
## m      No. of biomarker candidates                                                                                             ##
## k      No. of "true" biomarkers among the m biomarker candidates                                                               ##
## p.0    (Largest) p-value of the true biomarker candidates                                                                      ##
##                                                                                                                                ##
## Value                                                                                                                          ## 
## Probability as described above.                                                                                                ##                                   
####################################################################################################################################
prob.above.ci <- function(m,k,p.0,alpha=0.05){
  return(1 - pbinom(qbinom(1-alpha,m,p.0) - k,m-k,p.0))
}


####################################################################################################################################
## Computes the number of required "true" biomarkers with a p-value of p.0 (or smaller) among m biomarker candidates such that    ##
## the number of biomarker candidates exceeds the (1-alpha) CI with probability pow.                                              ##
##                                                                                                                                ##
## Arguments                                                                                                                      ##
## m      No. of biomarker candidates                                                                                             ##
## p.0    (Largest) p-value of the true biomarker candidates                                                                      ##
## alpha  Parameter to define the confidence level (1-alpha)                                                                      ##
## pow    Minimum probability with which the no of biomarker candidates the (1-alpha) CI                                          ##
##                                                                                                                                ##
## Value                                                                                                                          ## 
## Number of required "true" biomarkers with a p-value of p.0 (or smaller)                                                        ##    
####################################################################################################################################
min.no.biomarkers <- function(m,p.0,alpha=0.05,pow=0.8){
  k <- 0
  is.found <- F
  while(!is.found){
    k <- k+1
    if (prob.above.ci(m,k,p.0,alpha=alpha)>=pow){
      is.found <- T
    }
  }
  return(k)
}


####################################################################################################################################
## Plots the number of required "true" biomarkers depending on the number of biomarker candidates and the p-values.               ##
## (Compare min.no.biomarkers).                                                                                                   ##
##                                                                                                                                ##
## Arguments                                                                                                                      ##
## m      Array of numbers of biomarker candidates                                                                                ##
## p.0    Array of (largest) p-value of the true biomarker candidates                                                             ##
## alpha  As in min.no.biomarkers                                                                                                 ##
## pow    As in min.no.biomarkers                                                                                                 ##
##                                                                                                                                ##
## Value                                                                                                                          ## 
## Number of required "true" biomarkers with a p-value of p.0 (or smaller)                                                        ##
## xlab, ylab, zlab    Labels for the axes                                                                                        ##
## zlim   Range for the z-axis                                                                                                    ##
## theta, phi   Rotation angles for the 3D-plot                                                                                   ##  
## grid.plot  Indicates whether the surface should be endowed with a grid (TRUE) or with different colours (FALSE)                ##
## col.grid  Colour for the surface in case of a grid plot                                                                        ##  
####################################################################################################################################
plot.no.required.true.biomarkers <- function(m.values,p.values,alpha=0.05,pow=0.8,xlab="No. of biomarker candidates",ylab="p-value",zlab="Required no. of biomarkers",zlim=NULL,main=NULL,theta=30,phi=50,grid.plot=T,col.grid="cyan"){
  z <- matrix(0,nrow=length(m.values),ncol=length(p.values))
  for (i in 1:length(m.values)){
    for (j in 1:length(p.values)){
      z[i,j] <- min.no.biomarkers(m.values[i],p.values[j],alpha=alpha,pow=pow)
    }
  }
  # if (grid.plot){
  #   persp(m.values,p.values,z,theta=theta,phi=phi,r=2,shade=0.4,axes=T,scale=T,box=T,nticks=5,ticktype="detailed",col=col.grid,xlab=xlab,ylab=ylab,zlab=zlab,main=main)
  # }else{
  #   persp3D(m.values,p.values,z,theta=theta,phi=phi,axes=T,zlim=zlim,scale=2,box=TRUE,nticks=5,ticktype="detailed",xlab=xlab, ylab=ylab,zlab=zlab,main=main)
  # } 
  
  plot_ly(x= m.values,
          y= p.values,
          z= z,
          type = "contour" 
  )
}
