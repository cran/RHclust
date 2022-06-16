VIPnoSNP = function (Simulated = NULL, CPG = NULL, GE = NULL,
                CPGname = NULL, GEname = NULL, v,
                optimize = c('off','min','slope','elbow'),
                iter_max = 1000, nstart = 5, fit = c('aic','bic'),
                seed = NULL, ct = c('mean','median'), verbose = FALSE){
  if(length(optimize) > 1) optimize = 'off'
  if(length(fit) > 1) fit = 'aic'
  if(length(ct) > 1) ct = 'mean'

  #' Simulated data parse
  if(is.list(Simulated)){
    # Names
    cpg.name = as.data.frame(Simulated$CPG_Index);colnames(cpg.name) = c('CPG','gene')
    ge.name = as.data.frame(Simulated$GE_Index);colnames(ge.name) = c('GE','gene')

    obs.cpg = as.data.frame(Simulated$CPG);colnames(obs.cpg) = cpg.name[,1]
    obs.cpg = data.matrix(obs.cpg)

    obs.ge = as.data.frame(Simulated$GE);colnames(obs.ge) = ge.name[,1]
    obs.ge  = data.matrix(obs.ge)


    # Cluster Assignments
    SimClusAssign = data.matrix(unlist(Simulated[which(names(Simulated)=='Vec')]))
  }

  #' Read in data
  if(!is.null(CPGname)) {
    cpg.name = as.data.frame(CPGname);colnames(cpg.name) = c('CPG','gene')
  }
  if(!is.null(GEname))  {
    ge.name = as.data.frame(GEname);colnames(ge.name) = c('GE','gene')
  }


  if(!is.null(CPG)) {
    obs.cpg = as.data.frame(CPG);colnames(obs.cpg) = cpg.name[,1]
    obs.cpg = data.matrix(obs.cpg)
  }
  if(!is.null(GE))  {
    obs.ge = as.data.frame(GE);colnames(obs.ge) = ge.name[,1]
    obs.ge  = data.matrix(obs.ge)
  }



  #' Check Observed and Named data sets for NAs

  if(anyNA(obs.cpg) || anyNA(cpg.name)){
    cat('Missing values in Observed CPGs:',sum(anyNA(obs.cpg)),' \n')
    cat('Missing values in Named CPGs:',sum(anyNA(cpg.name)),' \n')
    stop('Either Observed or Named CPGs contain missing values \n')}

  if(anyNA(obs.ge) || anyNA(ge.name)){
    cat('Missing values in Observed GEs:',sum(anyNA(obs.ge)),' \n')
    cat('Missing values in Named GEs:',sum(anyNA(ge.name)),' \n')
    stop('Either Observed or Named GEs contain missing values \n')}

  #' Check Observed data sets for equals rows
  if (nrow(obs.ge) != nrow(obs.cpg))
    stop('Observed data must have the same number of rows \n')
  n = as.integer(nrow(obs.ge))

  #' Check cluster input
  if(missing(v))    stop('Number of clusters missing for v. \n')
  if(length(v) > 2) stop('v must be formated as either v = k or v = c(k1,k2). \n')
  if (max(v) > n)   stop('More clusters than number of rows. \n')
  v = as.integer(v); if (length(v) == 1) optimize = 'off' # No need to optimize if v == 1

  #' Check penalizer
  if (fit != 'aic' && fit != 'bic' && optimize != 'off') stop('Optimize must be "aic" or "bic". \n')

  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  elbow <- function(x, y) {

    max_df = data.frame(x = c(x[1], tail(x,1)), y = c(y[1], tail(y,1)))

    # Creating straight line between the max values
    fit <- lm(max_df$y ~ max_df$x)

    # Distance from point to line
    distances <- c()
    for(i in 1:length(x)) {
      distances <- c(distances,
                     abs(coef(fit)[2]*x[i] - y[i] + coef(fit)[1]) /
                       sqrt(coef(fit)[2]^2 + 1^2))
    }

    # Max distance point
    x_max_dist <- x[which.max(distances)]
    y_max_dist <- y[which.max(distances)]

    return(x_max_dist)
    #cat(c(x_max_dist, y_max_dist),'\n')
    #cat('Optimal clusters:',x_max_dist,'\n')
  }

  #' Identify upper and lower bounds if any
  if (length(v) != 1) {lb = v[1];ub=v[2]} else {lb = ub = v}
  if(lb>ub) stop('Lower bound must be smaller than the upper bound \n')

  # Initialize some outputs
  TOTaic =TOTbic = TOTwss = OPwss = OPaic = OPbic = cOPwss = output.all =slope.output = NULL
  #output.all=vector(mode = "list", length = ub)

  if (!is.null(seed)){
    set.seed(seed)
  }

  #' Start loop for range of clusters
  for (k in lb:ub){
    clus = k
    k_input = k


    #' Random sample indexes for k number of clusters
    sampid = sample(nrow(obs.ge), clus, replace = F)

    #' Rows of centers based on the indexes
    GEcenters = matrix(obs.ge[sampid, ], nrow = k); colnames(GEcenters) = colnames(obs.ge)
    CPGcenters = matrix(obs.cpg[sampid, ], nrow = k); colnames(CPGcenters) = colnames(obs.cpg)

    #' Used for nstarts

    #' Initalize variables for main loop
    clusters = n
    tot.dists = NULL
    moved = NULL
    iter = 1


    #' Main loop
    while (iter < iter_max){
      dists = matrix(0, nrow = n, ncol = k)


      for (g in unique(ge.name$gene)){
        #' Values only for the gth gene
        ge.index = (ge.name[which(ge.name$gene==g),])
        cpg.index = (cpg.name[which(cpg.name$gene==g),])

        for (i in 1:nrow(GEcenters)){
          #GE
          obs.ge1 = as.matrix(obs.ge[,g])
          gedist = (obs.ge1 - matrix(rep(GEcenters[i,g], nrow(obs.ge1)),
                                     nrow = nrow(obs.ge1), byrow= T))^2 # Euclidian Distance
          gedistsums = rowSums(gedist)


          #CPG
          obs.cpg1 = as.matrix(obs.cpg[,cpg.index[,1]])
          cpgdist = (obs.cpg1 - matrix(rep(CPGcenters[i,cpg.index[,1]], nrow(obs.cpg1)),
                                       nrow = nrow(obs.cpg1), byrow=T))^2 # Euclidian Distance
          cpgdistsum = rowSums(cpgdist)
          #cpgdistsum = cpgdist

          #Distance
          distance = sqrt(cpgdistsum + gedistsums)
          dists[,i] = dists[,i] + distance
        }
      }


      #cat(dim(GEcenters),dim(CPGcenters), dim(SNPcenters), '\n')
      #' Cluster assignments
      old.clusters = clusters
      if (k_input == 1){
        clusters = rep(1, n)
      } else {clusters = matrix(apply(dists[,colSums(dists != 0) > 0], 1, which.min))}
      moved = c(moved, sum(clusters != old.clusters))

      #' Update clusters by means for numerics and mode for categorical
      #'
      #' # Add option for FUN for mean or median
      GEcenters = as.matrix(aggregate(obs.ge, by=list(clusters), FUN=ct)[,-1]);colnames(GEcenters) = colnames(obs.ge)
      CPGcenters = as.matrix(aggregate(obs.cpg, by=list(clusters), FUN=ct)[,-1]);colnames(CPGcenters) = colnames(obs.cpg)


      #' If none of the values moved then we can break the loop
      if(moved[length(moved)] == 0) break

      #' Also the same as moved
      iter = iter+1

    }

    #' AIC and BIC penalizers and stores values
    #' Penalize based on the number of clusters

    #' Assign clusters based on min distance
    mins = apply(cbind(clusters, dists[,colSums(dists != 0) > 0]), 1, function(z) z[z[1] + 1])
    within = as.numeric(by(mins, clusters, sum))
    #if(length(unique(clusters))==1){within =as.numeric(aggregate(dists, by=list(clusters), FUN=sum)[,-1])}
    tot_within = sum(within)
    tot.dists = c(tot.dists, sum(tot_within))

    Centrow = k*2
    TOTwss[k] = tot_within
    TOTaic[k] = TOTwss[k] + (2*k*Centrow)
    TOTbic[k] = TOTwss[k] + (log(n)*Centrow)


    #' Normal output list
    output = list('Size' = table(clusters),
                  'Cluster' = clusters,
                  'GE_Centers' = t(GEcenters),
                  'CPG_Centers' = t(CPGcenters),
                  'WSS' = within,
                  'Total_WSS' = tot_within,
                  'Moved' = moved,
                  'AIC' = TOTaic[k],
                  'BIC' = TOTbic[k])

    if (is.list(Simulated)){
      BCtest = cbind(SimClusAssign,clusters)
      output[['BC Test']] = BCtest
    }

    #' Number of starts since KNN is relative to starting position.
    #' Highly recommend at least 5 nstarts
    if (nstart > 1){
      for (j in 2:nstart){
        if (!is.null(seed)){seed=j}
        output.new = VIPnoSNP(Simulated = Simulated, CPG = CPG, GE = GE,
                         CPGname = CPGname, GEname = GEname,
                         v = k_input, iter_max = iter_max, nstart = 1,
                         seed = seed, ct = ct)
        if(((output.new$Total_WSS < output$Total_WSS))&&(output.new$Total_WSS !=0)) output = output.new
        if(j == nstart && (length(output$Size) < k_input)){
          warning(k_input,' initial clusters reduced to ',(length(output$Size)),'\n')
        }
      }
    }

    if(optimize != 'off'){output.all[[k]]=output}
    cOPwss[k] = length(output$Size)
    OPwss[k] = output$Total_WSS
    OPaic[k] = output$AIC;
    OPbic[k] = output$BIC;

    if (optimize == 'slope'){
      PenType = switch(ifelse(fit == 'aic',1,2),
                       OPaic,
                       OPbic)
      if(sum(!is.na(PenType)) > 2){
        #' If the difference in point b and a are > 0 then the slope is increasing
        if(((PenType[k] - PenType[k-1]) > 0)){ #&& is.null(slope.output)){
          #' Cluster is assigning based on the first instance of an increasing slope
          #'
          if (verbose){
            cat('Optimize Type:',optimize, '\n')
            cat('    Penalizer:',fit,'\n')
            cat('     Clusters:',k-1,' \n')
          }

          opt.output = output.all[[k-1]]
          return(opt.output)


        }
      }
    }
  } #' End of loop for l : u


  #' Used for assigning which penalizing factor to use for optimal cluster assignment
  PenType = switch(ifelse(fit == 'aic',1,2),
                   OPaic,
                   OPbic)

  if(optimize == 'min') {
    opt.output=output.all[[which.min(PenType)]]
  }

  if(optimize == 'elbow'){
    elb = elbow(lb:ub,PenType[!is.na(PenType)])
    opt.output = output.all[[elb]]
  }


  #' Creates plot for the Total WSS, AIC, and BIC to check for user defined optimal clusters
  if (length(v) != 1){
    outputPlot = list('Total_WSS' = OPwss,
                      'AIC' = OPaic,
                      'BIC' = OPbic,
                      'Clusters' = cOPwss)

    #' Partition all 3 plots into one frame
    par_ = par(no.readonly = TRUE)
    on.exit(par(par_))
    par(mfrow = c(1,3))

    outputPlot = lapply(outputPlot, na.omit)

    plot(x = lb:ub,outputPlot$Total_WSS, type = 'b',
         ylim = c(min(outputPlot$Total_WSS), max(outputPlot$Total_WSS)*1.1),
         xlab = 'Clusters', ylab = '', main = 'Total_WSS')
    if (any(seq(lb:ub) != outputPlot$Clusters)){
      text(x = lb:ub, y = outputPlot$Total_WSS, labels = paste0('(',outputPlot$Clusters,')'),pos = 3, col = 'red')
      warning('Total WSS clusters have been labeled based on assigned number of centers.')
    }

    plot(x = lb:ub, outputPlot$AIC, type = 'b',
         ylim = c(min(outputPlot$AIC), max(outputPlot$AIC)*1.1),
         xlab = 'Clusters', ylab = '', main = 'AIC')
    if (any(seq(lb:ub) != outputPlot$Clusters)){
      text(x = lb:ub,outputPlot$AIC,labels = paste0('(',outputPlot$Clusters,')'),pos = 3, col = 'red')
      warning('AIC clusters have been labeled based on assigned number of centers.')
    }

    plot(x = lb:ub, outputPlot$BIC, type = 'b',
         ylim = c(min(outputPlot$BIC), max(outputPlot$BIC)*1.1),
         xlab = 'Clusters', ylab = '', main = 'BIC')
    if (any(seq(lb:ub) != outputPlot$Clusters)){
      text(x = lb:ub,outputPlot$BIC,labels = paste0('(',outputPlot$Clusters,')'),pos = 3, col = 'red')
      warning('BIC clusters have been labeled based on assigned number of centers.')
    }

    output[["PlotData"]]=outputPlot
    if(optimize == 'off') return(outputPlot)
  }
  if(optimize != 'off'){
    if(optimize == 'slope') {
      warning('No increasing slope found. Returning minimum fit criteria and largest cluster. \n')
      opt.output = output.all[[which.min(PenType)]]
    }

    if (verbose){
      cat('Optimize Type:',optimize, '\n')
      cat('    Penalizer:',fit,'\n')
      cat('     Clusters:',length(opt.output$Size),' \n')
    }

    opt.output[["PlotData"]]=outputPlot
    return(opt.output)
  }
  return(output) # Simply returns the normal output list plus plots if more than one value for cluster provided

}

