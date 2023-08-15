require(Runuran)

SimData = function(seed = NULL, gene = 36,
                   k = c(33,33,34),
                   GEbar = 5, GEsd = 0.5,
                   CPGbar = 4, CPGsd = 0.5, SameCPG = FALSE,
                   SameSNP = FALSE, ProbDist = NULL, SameGeneDist = TRUE, distinct=TRUE){
  #require(Runuran)
  GE = CPG = SNP = gsd = SNP_idx = CPG_idx = Vec = NULL

  # AGE = BMI = SYSBP = DYSBP = GROUP = GROUPCAT = NULL
  # AGEm = BMIm = SYSBPm = DYSBPm = 0
  #
  # AGEbar = 30
  # BMIbar = 18.5
  # SYSBPbar = 120
  # DYSBPbar = 80


  gem = cpgm = 0
  a = b = floor(gene/3); (c = gene - a - b)
  if(gene < 3) a = b = c = 1
  o = list(1:a,
           a+1:b,
           a+b+1:c)
  problist = ProbDist
  if (is.null(problist)){
    problist = list(c(0.50,0.25,0.25),
                    c(0.20,0.55,0.25),
                    c(0.30,0.15,0.55),
                    c(0.20,0.50,0.30),
                    c(0.45,0.20,0.35),
                    c(0.26,0.52,0.22),
                    c(0.23,0.21,0.56),
                    c(0.23,0.52,0.25),
                    c(0.51,0.25,0.24),
                    c(0.21,0.58,0.21))
  }
  if(length(k) > length(problist)) stop('Number of clusters exceeds ProbDist.\n')
  probdist = 1

  for (g in 1:gene){
    set.seed(g)
    sd = round(runif(1,3,50))
    gsd = rbind(gsd,sd)
    SNP_idx = rbind(SNP_idx, do.call('cbind',list(replicate(sd,paste0('GE',g)))))
  }
  set.seed(seed)
  for (i in 1:length(k)){
    Vec = c(Vec, rep(i,k[i]))


    GE = rbind(GE,replicate(n = gene, Runuran::urnorm(k[i], GEbar+gem, sd = GEsd)))

    CPG = rbind(CPG, t(replicate(k[i],Runuran::urnorm(n = (15*(a+c)), mean = (CPGbar+cpgm), sd = CPGsd))))

    SNP = rbind(SNP,t(replicate(k[i],sample(LETTERS[1:3],sum(gsd), replace=TRUE, prob=problist[[probdist]]))))

    cpgm = cpgm + 4
    gem = gem + 5
    probdist = i
    if(SameCPG) cpgm = 0
    if(SameSNP) probdist = 1

  }

  for (cpg_idx in 1:3){
    CPG_idx = rbind(CPG_idx,data.frame(rep(paste0('GE',o[[cpg_idx]]), cpg_idx*5)))
  }

  # if (isTRUE(cov)){
  #   for (i in 1:length(k)){
  #     # Numerics
  #     AGE = rbind(AGE,replicate(n = 1, Runuran::urnorm(k[i], AGEbar+AGEm, sd = 1)))
  #     AGEbar = AGEbar + 7
  #
  #     BMI = rbind(BMI,replicate(n = 1, Runuran::urnorm(k[i], BMIbar+BMIm, sd = 1)))
  #     BMIm = BMIm + 7
  #
  #     SYSBP = rbind(SYSBP,replicate(n = 1, Runuran::urnorm(k[i], SYSBPbar+SYSBPm, sd = 1)))
  #     SYSBPm = SYSBPm + 14
  #
  #     DYSBP = rbind(DYSBP,replicate(n = 1, Runuran::urnorm(k[i], DYSBPbar+DYSBPm, sd = 1)))
  #     DYSBPm = DYSBPm + 8
  #
  #     # Categorical
  #     GROUP = c(GROUP,matrix(rep(i,k[i])))
  #
  #     GROUPCAT = c(GROUPCAT, matrix(rep(LETTERS[i],k[i])))
  #   }
  #
  #   GROUP = as.factor(GROUP)
  #   GROUPCAT = as.factor(GROUPCAT)
  #   #cov = data.frame(AGE,BMI,SYSBP,DYSBP,GROUP,GROUPCAT)
  # }
  #
  # COV = data.frame(AGE,BMI,SYSBP,DYSBP,GROUP,GROUPCAT)
  tbl = NULL
  NMS = NULL
  for(cv in 1:10) {
    nums <- paste("NUM", cv, sep = "")
    cats <- paste("CAT", cv, sep = "")

    means = seq(10,100,10)
    val = grp = 0
    meandelt = 0
    covsd = 1
    if(SameGeneDist){
      for (j in 1:length(k)){
        #val = rbind(val,replicate(n = 1, Runuran::urnorm(k[j], means[cv]+meandelt, sd = covsd)))
        val = rbind(val,replicate(n = 1, runif(k[j], means[cv]+meandelt -1 , means[cv]+meandelt + 1)))
        #meandelt = meandelt + (means[cv]/2)

        #cat('Means: ',means[cv], 'MeanDeltLower: ',means[cv]+meandelt -1, 'MeanDeltUpper: ',means[cv]+meandelt + 1,'\n')

        grp = c(grp,matrix(rep(means[cv]+floor(meandelt),k[j])))
        meandelt = meandelt + (means[cv]/2)

      }
    }
    if(!SameGeneDist){
      for (j in 1:length(k)){
        val = rbind(val,replicate(n = 1, Runuran::urnorm(k[j], means[cv]+meandelt, sd = covsd)))
        #meandelt = meandelt + (means[cv]/2)

        #cat('Means: ',means[cv]+meandelt, 'SD: ', covsd,'\n')

        grp = c(grp,matrix(rep(means[cv]+floor(meandelt),k[j])))
        meandelt = 0 #meandelt + (means[cv]/2)


      }
    }

    meandelt = 0


    #grp = c(grp,matrix(rep(cv,k[j])))


    tbl = cbind(tbl,assign(nums,val[-1]))

    assign(cats,(grp[-1]))

    tbl = cbind(tbl,assign(cats,(grp[-1])))
    NMS = c(NMS,nums,cats)
  }


  tbl = as.data.frame(tbl)
  catcols = c(2,4,6,8,10,12,14,16,18,20)
  tbl[catcols] = lapply(tbl[catcols],factor)
  colnames(tbl) = NMS
  COV = tbl[,order(names(tbl))]

  COV = COV[,c(1,3,4,5,6,7,8,9,10,2,11,13,14,15,16,17,18,19,20,12)]


  colnames(GE) = paste0("GE",1:ncol(GE)) ;
  colnames(CPG) = paste0("CPG",1:ncol(CPG));
  colnames(SNP) = paste0("SNP",1:ncol(SNP))

  GEn = cbind(paste0("GE",1:ncol(GE)),paste0("GE",1:ncol(GE)));
  colnames(GEn) = c('GE','gene')

  CPGn = cbind(paste0('CPG',1:ncol(CPG)), CPG_idx);
  colnames(CPGn) = c('CPG','gene')

  SNPn = cbind(paste0("SNP",1:ncol(SNP)),SNP_idx);
  colnames(SNPn) = c('SNP','gene')

  if(gene == 1) {CPGn = CPGn[1:5,]; CPG = CPG[,1:5]}
  if(gene == 2) {CPGn = CPGn[1:15,]; CPG = CPG[,1:15]}


  res = list('Clusters' = k,
             'Vec' = sort(Vec),
             'GE' = data.frame(GE),
             'CPG' = data.frame(CPG),
             'SNP' = data.frame(SNP),
             'GE_Index' = data.frame(GEn),
             'CPG_Index' = data.frame(CPGn),
             'SNP_Index' = data.frame(SNPn),
             'Covariates' = data.frame(COV))

  if (distinct){
    NUM1 = matrix(res$Covariates[[11]])
    NUM2 = NUM1*3
    NUM3 = NUM2*3
    NUM4 = NUM3*3
    NUM5 = NUM4*3
    NUM6 = NUM5*3
    NUM7 = NUM6*3
    NUM8 = NUM7*3
    NUM9 = NUM8*3
    NUM10 = NUM9*3
    newcovs = cbind(NUM1,NUM2, NUM3,NUM4,NUM5,NUM6,NUM7,NUM8,NUM9,NUM10)/10
    res$Covariates = data.frame(cbind(res$Covariates[1:10],newcovs))
  }

  if (!SameGeneDist & distinct){
    cats = replicate(10,sample(1:length(k),sum(k), replace=TRUE))
    means = c(10,40,80,120,160,200,240,280,320,360)
    val = grp = 0
    meandelt = 0
    covsd = 5
    #k = c(33,33,34)
    nums = NULL
    val = NULL
    for(cv in 1:10) {
      #cv = 1
      for (j in 1:length(k)){
        val = rbind(val,replicate(n = 1, Runuran::urnorm(k[j], means[cv]+meandelt, sd = covsd)))
        #meandelt = meandelt + (means[cv]/2)

        #cat('Means: ',means[cv]+meandelt, 'SD: ', covsd,'\n')

        grp = c(grp,matrix(rep(means[cv]+floor(meandelt),k[j])))
        meandelt = 0 #meandelt + (means[cv]/2)

      }
      nums = cbind(nums,val)
      #nm = cbind(nm,val)
      val = NULL
    }

    #nums = res$Covariates[c(11:20)]
    res$Covariates[c(11:20)] = nums

    res$Covariates = as.data.frame(cbind(cats,nums))
    res$Covariates[c(1:10)] = lapply(res$Covariates[c(1:10)],factor)
  }
  if (!SameGeneDist & !distinct){
    cats = replicate(10,sample(1:length(k),sum(k), replace=TRUE))
    means = seq(10,100,10)
    val = grp = 0
    meandelt = 0
    covsd = 5
    #k = c(33,33,34)
    nums = NULL
    val = NULL
    for(cv in 1:10) {
      #cv = 1
      for (j in 1:length(k)){
        val = rbind(val,replicate(n = 1, Runuran::urnorm(k[j], means[cv]+meandelt, sd = covsd)))
        #meandelt = meandelt + (means[cv]/2)

        #cat('Means: ',means[cv]+meandelt, 'SD: ', covsd,'\n')

        grp = c(grp,matrix(rep(means[cv]+floor(meandelt),k[j])))
        meandelt = 0 #meandelt + (means[cv]/2)

      }
      nums = cbind(nums,val)
      #nm = cbind(nm,val)
      val = NULL
    }

    #nums = res$Covariates[c(11:20)]
    res$Covariates[c(11:20)] = nums

    res$Covariates = as.data.frame(cbind(cats,nums))
    res$Covariates[c(1:10)] = lapply(res$Covariates[c(1:10)],factor)
  }

  return(res)
}


