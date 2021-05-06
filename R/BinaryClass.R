BinaryClass = function(x){

  if(class(x)[1] != 'table' && ncol(x) > 2) stop('x must be a table or data.frame of 2 columns')
  if(class(x)[1] != 'table'){tab = table(x[,2], x[,1])}

  if(class(x)[1] == 'table'){tab = x}


  # Important stuff
  acc = NULL
  if (nrow(tab) == 1){
    m = tab[which.max(tab)]
    t = which.max(tab)
    n = sum(tab)
    acc = m/n
  }
  m = apply(tab,2,max)
  t = apply(tab, 2, which.max)
  cs = colSums(tab)
  rs = rowSums(tab)
  n = sum(tab)

  subs = matrix(0, nrow = nrow(tab), ncol = ncol(tab))
  dupmodes = as.numeric(t[duplicated(t)])

  if (any(duplicated(t))){

    if (length(unique(dupmodes)) == 1){
      subs = sum(tab[unique(dupmodes),]) - m[which.max(tab[dupmodes,])]
    } else subs = rowSums(tab[dupmodes,]) - m[which.max(tab[,dupmodes])]

  }

  #if (all(subs == 0) && !is.na(subs)){subs = 0}
  if (sum(subs) == 0 || is.na(sum(subs))){subs = 0}

  if (is.null(acc)){acc = (sum(m)-sum(subs)) / n}

  SE = ((m) / cs)

  # Specificity
  SP = (n + (m) - cs - rs[t]) / (n - cs)

  # Prevelance
  Prev = cs / n

  # PPV
  PPV = (SE * Prev) / ((SE*Prev) + ((1 - SP) * (1 - Prev)))

  # NPV
  NPV = (SP * (1-Prev)) / (((1 - SE) * Prev) + ((SP)*(1-Prev)))

  # Detection Rate
  DR = m / n

  # Detection Prev
  DP = rs[t] / n

  # Balanced Accuracy
  BA = (SE + SP)/2


  # Clopper-Pearson Interval for Confidence Interval
  x = acc * n
  lb = (1+ (n - x + 1)/ (x*qf(.025,2*x, 2*(n - x + 1))))^-1

  if(n - x == 0){ub = 1}
  else ub = (1 + (n - x) / ((x + 1)* qf(.975,2*(x+1), 2*(n-x))))^-1
  CI = c(lb,ub)

  # Outputs
  names(dimnames(tab)) = c('predicted', 'true')
  out = round(rbind(SE, SP, PPV, NPV, Prev, DR, DP, BA),4)
  rownames(out) = c('Sensitivity','Specificity', 'Pos Pred Value', 'Neg Pred Value',
                    'Prevalence', 'Detection Rate', 'Detection Prev', 'Balanced Accuracy')
  colnames(out) = paste0('Class ',seq(1:ncol(tab)))
  res = list('Table' = tab,
             'Accuracy' = acc,
             '95% CI' = CI,
             'Group Measures' = out)
  return(res)
}

