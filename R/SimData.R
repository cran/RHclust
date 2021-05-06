require(Runuran)

SimData = function(seed = NULL, gene = 36,
                   k = c(33,33,34),
                   GEbar = 5, GEsd = 0.5,
                   CPGbar = 4, CPGsd = 0.5, SameCPG = FALSE,
                   SameSNP = FALSE, ProbDist = NULL){
  #require(Runuran)
  GE = CPG = SNP = gsd = SNP_idx = CPG_idx = Vec = NULL
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
                    c(0.45,0.20,0.35))
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
             'SNP_Index' = data.frame(SNPn))
  return(res)
}

#' SimData = function(s = NULL, g = 36, k = c(33,33,34), type = 'Default'){
#'
#'
#'   # Required Packages
#'   require(Runuran)
#'   require(plyr)
#'   require(data.table)
#'   #' Require at least 3 genes because we divide it by 3
#'   #if (g < 3) stop('g must be 3 or greater. \n')
#'   if(sum(k) != 100) stop('Sum of clusters must equal 100')
#'   og = g
#'   if (og < 3){
#'     g = 3
#'   }
#'   gene = g
#'
#'   #' Number of genes per cluster
#'   clusters = k
#'
#'   #' Why we need at least 3 genes
#'   #' This just breaks up the genes for the CPG assignments
#'   if (g == 1 || g == 2){
#'     a_ = b_ = c_ = 1
#'     cpgs = c(a_,b_,c_)
#'   } else{a_ = b_ = floor(gene/3); (c_ = gene - a_ - b_)
#'   cpgs = c(a_,b_,c_)}
#'   # a_ = b_ = floor(gene/3); (c_ = gene - a_ - b_)
#'   # cpgs = c(a_,b_,c_)
#'
#'     # Probability Distribution for SNPs
#'     if (type == 'SameSNP'){
#'       prob.c.1<-c(0.50,0.25,0.25)
#'       prob.c.2<-c(0.50,0.25,0.25)
#'       prob.c.3<-c(0.50,0.25,0.25)
#'       prob.c.4<-c(0.50,0.25,0.25)
#'     } else{
#'       prob.c.1<-c(0.50,0.25,0.25)
#'       prob.c.2<-c(0.20,0.55,0.25)
#'       prob.c.3<-c(0.30,0.15,0.55)
#'       prob.c.4<-c(0.20,0.50,0.30)
#'     }
#'
#'
#'   #' Initialize some values
#'   GE = CPG = CPG_n = SNPprev = SNP = SNPs = Vec = NULL
#'   b = 1
#'
#'   m = 1:a_#(1:(gene/length(clusters)))
#'
#'   # if (!is.null(s)){
#'   #   set.seed(s)
#'   # }
#'
#'   for (j in clusters){
#'
#'     clusid = b
#'
#'     # SNP
#'     gsd = NULL
#'     SNP_n = NULL
#'     for (g in 1:gene){
#'       set.seed(g)
#'       sd = round(runif(1,3,50))
#'       gsd = rbind(gsd,sd)
#'       SNP_n = rbind(SNP_n,do.call('cbind',list(replicate(sd,paste0('GE',g)))))
#'     }
#'
#'     set.seed(s)
#'
#'     #set.seed(s)
#'     SNPs = matrix(sample(LETTERS[1:3], sum(gsd)*j, replace=TRUE, prob=get(paste0('prob.c.',(clusid)))),nrow=j,ncol= sum(gsd))
#'     SNP = rbind(SNP,SNPs)
#'
#'     # GE
#'     GEs = matrix(urnorm(n = gene*j, mean = (clusid*5), sd=1, lb = 0, ub = Inf), nrow = j, ncol = gene)
#'     GE = rbind(GE,GEs)
#'
#'     #' CPG
#'     #' a_ = b_
#'     #' n = 5a_ + 10b_ + 15c_ ---> n = (5a_ + 10a_) + 15c_ ---> 15a_ + 15c_
#'     #' n = 15(a_ + c_)
#'     if (type == 'SameCPG'){
#'       CPGs = matrix(urnorm(n = (15*(a_+c_))*j, mean = 4, sd = 0.5),
#'                     nrow = j, ncol = 15*(a_+c_))
#'     } else {CPGs = matrix(urnorm(n = (15*(a_+c_))*j, mean = (clusid*4), sd = 0.5),
#'                           nrow = j, ncol = 15*(a_+c_))}
#'
#'     CPG = rbind(CPG,CPGs)
#'
#'     CPG_n = rbind(CPG_n,data.frame(rep(paste0('GE',m), clusid*5)))
#'
#'     Vec = rbind(matrix(rep(b, k[b])), Vec)
#'
#'     if (b==1)m = m + b_#m + (gene/length(clusters))
#'     if (b==2)m = (tail(m, n = 1)+1):gene #(m + b_)
#'     b = b + 1
#'
#'   }
#'   colnames(GE) = paste0("GE",1:ncol(GE)) ;
#'   colnames(CPG) = paste0("CPG",1:ncol(CPG));
#'   colnames(SNP) = paste0("SNP",1:ncol(SNP))
#'
#'   GEn = cbind(paste0("GE",1:ncol(GE)),paste0("GE",1:ncol(GE)));
#'   colnames(GEn) = c('GE','gene')
#'
#'   CPGn = cbind(paste0('CPG',1:ncol(CPG)), CPG_n);
#'   colnames(CPGn) = c('CPG','gene')
#'
#'   SNPn = cbind(paste0("SNP",1:ncol(SNP)),SNP_n);
#'   colnames(SNPn) = c('SNP','gene')
#'
#'
#'   res = list('Clusters' = k,
#'              'Vec' = sort(Vec),
#'              'GE' = data.frame(GE),
#'              'CPG' = data.frame(CPG),
#'              'SNP' = data.frame(SNP),
#'              'GE_Index' = data.frame(GEn),
#'              'CPG_Index' = data.frame(CPGn),
#'              'SNP_Index' = data.frame(SNPn))
#'
#'   if (og >= 3){
#'     return(res)
#'   }
#'
#'   if(og < 3){
#'
#'     GEn = res$GE_Index[1,]
#'     CPGn = (res$CPG_Index[which(res$CPG_Index$gene=='GE1'),])
#'     SNPn = res$SNP_Index[which(res$SNP_Index$gene == 'GE1'),]
#'
#'     GE = res$GE[,GEn$GE]
#'     CPG = res$CPG[,CPGn$CPG]
#'     SNP = res$SNP[,SNPn$SNP]
#'
#'     if (og == 2){
#'       GEn_2 = res$GE_Index[2,]
#'       CPGn_2 = (res$CPG_Index[which(res$CPG_Index$gene=='GE2'),])
#'       SNPn_2 = res$SNP_Index[which(res$SNP_Index$gene == 'GE2'),]
#'
#'       GE_2 = res$GE[,GEn_2$GE]
#'       CPG_2 = res$CPG[,CPGn_2$CPG]
#'       SNP_2 = res$SNP[,SNPn_2$SNP]
#'
#'       GEn = rbind(GEn, GEn_2)
#'       CPGn = rbind(CPGn, CPGn_2)
#'       SNPn = rbind(SNPn, SNPn_2)
#'
#'       GE = cbind(GE, GE_2)
#'       CPG = cbind(CPG, CPG_2)
#'       SNP = cbind(SNP, SNP_2)
#'
#'     }
#'
#'     if ((og == 1) & length(CPGn != 5)) CPGn = CPGn[1:5,]
#'
#'     res = list('Clusters' = k,
#'                'Vec' = sort(Vec),
#'                'GE' = data.frame(GE),
#'                'CPG' = data.frame(CPG),
#'                'SNP' = data.frame(SNP),
#'                'GE_Index' = data.frame(GEn),
#'                'CPG_Index' = data.frame(CPGn),
#'                'SNP_Index' = data.frame(SNPn))
#'   }
#' }
#'
