BinaryClass = function(x){

  n = x
  #original=c(rep(1,33),rep(2,33),rep(3,34))
  original=n[,1]

  tab = table(n[,2],n[,1])
  total = sum(tab)
  final.correct.cluster<-NULL
  cluster.all<-NULL
  OG = original
  tbl = NULL
  maxes = NULL
  maxidx = NULL
  SE = NULL
  SP = NULL
  PREV = NULL
  PPV = NULL
  NPV = NULL
  DETR = NULL
  DETP = NULL
  BA = NULL
  cluster.temp<-rep(0,nrow(n))
  study.temp<-rep(0,nrow(n))
  dup<-NULL
  tbl_ = 1

  for (tb in unique(n[,1])){
    ogtbl = table(n[n[,1] == tb,2])

    #tbl = rbind(tbl,ogtbl)


    if(length(colnames(t(ogtbl))[which(ogtbl == max(ogtbl))])==1){
      maxidx[tbl_] = colnames(t(ogtbl))[which(ogtbl == max(ogtbl))]
      maxes[tbl_] = max(ogtbl)

    }else{
      #maxidx[tbl_] = sample(colnames(t(ogtbl))[which(ogtbl == max(ogtbl))],1)
      maxidx[tbl_] = colnames(t(ogtbl))[which(ogtbl == max(ogtbl))][1]

    }


    tbl_ = tbl_ + 1
  }

  center = maxidx
  max.all = maxes

  if(any(duplicated(center)==TRUE)){dup<-which.max(max.all[duplicated(center)])}
  #if(any(duplicated(center)==TRUE) & length(dup) > 1){dup<-sample(dup,1)}
  if(any(duplicated(center)==TRUE) & length(dup) > 1){dup<-dup[1]}

  for (tb in unique(n[,1])){
    cluster.temp[which(n[,2] == center[tb])]<-tb
  }


  if(any(duplicated(center)==TRUE)){cluster.temp[which(n[,2] == center[dup])]<-rep(dup,length(which(n[,2] == center[dup])))}

  cluster<-matrix(rep(0,nrow(n)*length(unique(OG))),nrow(n),length(unique(OG)))


  for(i in unique(n[,1])){
    for(ii in 1:nrow(n)){
      if(OG[ii] == i && cluster.temp[ii] == i){cluster[ii,i]<-'TP'}
      if(OG[ii] == i && cluster.temp[ii] != i){cluster[ii,i]<-'FN'}
      if(OG[ii] != i && cluster.temp[ii] == i){cluster[ii,i]<-'FP'}
      if(OG[ii] != i && cluster.temp[ii] != i){cluster[ii,i]<-'TN'}
    }

    ses = (length(cluster[which(cluster[,i] == 'TP') ,i]))/
      (length(cluster[which(cluster[,i] == 'TP') ,i])+length(cluster[which(cluster[,i] == 'FN') ,i]))
    SE<-c(SE,ses)

    sps = length(cluster[which(cluster[,i] == 'TN') ,i])/
      (length(cluster[which(cluster[,i] == 'TN') ,i])+length(cluster[which(cluster[,i] == 'FP') ,i]))
    SP<-c(SP,sps)

    prevs = (length(cluster[which(cluster[,i] == 'TP') ,i]) + length(cluster[which(cluster[,i] == 'FN') ,i]))/
      (length(cluster[which(cluster[,i] == 'TP') ,i]) + length(cluster[which(cluster[,i] == 'TN') ,i]) +
         length(cluster[which(cluster[,i] == 'FP') ,i]) + length(cluster[which(cluster[,i] == 'FN') ,i]))

    PREV = c(PREV, prevs)

    #ppvs = (ses*prevs)/((ses*prevs) + ((1-sps)*(1-prevs)))
    ppvs = length(cluster[which(cluster[,i] == 'TP') ,i]) /
      (length(cluster[which(cluster[,i] == 'TP') ,i]) + length(cluster[which(cluster[,i] == 'FP') ,i]))
    PPV = c(PPV, ppvs)

    #npvs = (sps*(1-prevs))/(((1-ses)*prevs) + ((sps)*(1-prevs)))
    npvs = length(cluster[which(cluster[,i] == 'TN') ,i]) /
      (length(cluster[which(cluster[,i] == 'TN') ,i]) + length(cluster[which(cluster[,i] == 'FN') ,i]))
    NPV = c(NPV,npvs)

    detr = length(cluster[which(cluster[,i] == 'TP') ,i]) /
      (length(cluster[which(cluster[,i] == 'TP') ,i]) + length(cluster[which(cluster[,i] == 'TN') ,i]) +
         length(cluster[which(cluster[,i] == 'FP') ,i]) + length(cluster[which(cluster[,i] == 'FN') ,i]))
    DETR = c(DETR,detr)

    detp = (length(cluster[which(cluster[,i] == 'TP') ,i]) + length(cluster[which(cluster[,i] == 'FP') ,i])) /
      (length(cluster[which(cluster[,i] == 'TP') ,i]) + length(cluster[which(cluster[,i] == 'TN') ,i]) +
         length(cluster[which(cluster[,i] == 'FP') ,i]) + length(cluster[which(cluster[,i] == 'FN') ,i]))
    DETP = c(DETP,detp)

    bas = (ses + sps)/2
    BA = c(BA,bas)

  }


  # for (SEs in unique(n[,1])){
  #   ses = length(cluster[which(cluster[,SEs] == 'TP') ,SEs])/
  #     (length(cluster[which(cluster[,SEs] == 'TP') ,SEs])+length(cluster[which(cluster[,SEs] == 'FN') ,SEs]))
  #   SE<-c(SE,ses)
  #
  #   sps = length(cluster[which(cluster[,SEs] == 'TN') ,SEs])/
  #     (length(cluster[which(cluster[,SEs] == 'TN') ,SEs])+length(cluster[which(cluster[,SEs] == 'FP') ,SEs]))
  #   SP<-c(SP,sps)
  #
  #   prevs = (length(cluster[which(cluster[,SEs] == 'TP') ,SEs]) + length(cluster[which(cluster[,SEs] == 'FN') ,SEs]))/
  #     (length(cluster[which(cluster[,SEs] == 'TP') ,SEs]) + length(cluster[which(cluster[,SEs] == 'TN') ,SEs]) +
  #     length(cluster[which(cluster[,SEs] == 'FP') ,SEs]) + length(cluster[which(cluster[,SEs] == 'FN') ,SEs]))
  #
  #   PREV = c(PREV, prevs)
  #
  #   ppvs = (ses*prevs)/((ses*prevs) + ((1-sps)*(1-prevs)))
  #   PPV = c(PPV, ppvs)
  #
  #   npvs = (sps*(1-prevs))/(((1-ses)*prevs) + ((sps)*(1-prevs)))
  #   NPV = c(NPV,npvs)
  #
  #   detr = length(cluster[which(cluster[,SEs] == 'TP') ,SEs]) /
  #     (length(cluster[which(cluster[,SEs] == 'TP') ,SEs]) + length(cluster[which(cluster[,SEs] == 'TN') ,SEs]) +
  #     length(cluster[which(cluster[,SEs] == 'FP') ,SEs]) + length(cluster[which(cluster[,SEs] == 'FN') ,SEs]))
  #   DETR = c(DETR,detr)
  #
  #   detp = (length(cluster[which(cluster[,SEs] == 'TP') ,SEs]) + length(cluster[which(cluster[,SEs] == 'FP') ,SEs])) /
  #     (length(cluster[which(cluster[,SEs] == 'TP') ,SEs]) + length(cluster[which(cluster[,SEs] == 'TN') ,SEs]) +
  #     length(cluster[which(cluster[,SEs] == 'FP') ,SEs]) + length(cluster[which(cluster[,SEs] == 'FN') ,SEs]))
  #   DETP = c(DETP,detp)
  #
  #   bas = (ses + sps)/2
  #   BA = c(BA,bas)
  # }

  correctcluster<-rep(0,nrow(n))
  for(i in 1:nrow(n)){
    if(original[i] == cluster.temp[i]){correctcluster[i]<-1}
  }

  final.correct.cluster<-c(final.correct.cluster,sum(correctcluster)/nrow(n))
  cluster.all<-c(cluster.all,length(unique(n[,2])))

  ACC = final.correct.cluster

  # Clopper-Pearson Interval for Confidence Interval
  x = ACC * total
  lb = (1+ (total - x + 1)/ (x*qf(.025,2*x, 2*(total - x + 1))))^-1

  if(total - x == 0){
    ub = 1
  } else {ub = (1 + (total - x) / ((x + 1)* qf(.975,2*(x+1), 2*(total-x))))^-1}
  CI = c(lb,ub)

  names(dimnames(tab)) = c('predicted', 'true')
  out = list('Table' = tab,
             'Accuracy' = ACC,
             '95% CI' = CI,
             'Sensitivity' = SE,
             'Specificity' = SP,
             'Prevalence' = PREV,
             'PPV' = PPV,
             'NPV' = NPV,
             'Detection Rate' = DETR,
             'Detection Prev' = DETP,
             'Balanced Accuracy' = BA)

  return(out)
}
