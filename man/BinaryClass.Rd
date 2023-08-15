\name{BinaryClass}
\alias{BinaryClass}
\title{Binary Classification}
\description{
A confusion matrix but allows for anaylsis of non-equal level data classifications.
}
\usage{
BinaryClass(x)
}
\arguments{
\item{x}{Can be a data frame dimensions at least 2 rows and 2 columns meant to represent observed and predicted values where the observed (true) values are in the first column and predicted columns in the second column.}
}

\details{
BinaryClass() is similar to a confusion matrix with binary classification outputs. The true positive values per column are identified based on the maximum number of assignments per category.
}

\value{
\item{Table}{the results of 'table()' on 'x'}

\item{Accuracy}{overall accuracy of classification}

\item{CI}{confidence interval of overall accuracy using Clopper-Pearson Interval}

\item{Group Measures}{the sensitivity, specificity, positive predictive value, negative predictive value, prevelance detection rate, detection prevalence, and balanced accuracy for each class}
}

\author{
jkhndwrk@memphis.edu
}

\examples{

# Basic example
true = c(rep(1,5), rep(2,5), rep(3,5), rep(4,5))
pred = c(rep(1,4),4,rep(2,5),2,rep(3,4),1,rep(4,4))
df = cbind(true,pred)
BinaryClass(df)


true = c(rep(1,5), rep(2,5), rep(3,5), rep(4,5))
pred = c(rep(1,5),rep(2,5),rep(3,10))
df = cbind(true,pred)
BinaryClass(df)


\donttest{
sd = SimData(k = c(10,40,50))
out = VIP(sd, v = 3, optimize = 'elbow', nstart = 5)
df = out$`BC Test`
BinaryClass(df)


## Looping through different clusters

sd = SimData(seed = 1, gene = 1)
acc = NULL
for (i in 1:5){
 out = VIP(sd, v = i, optimize = 'off', nstart = 5)
 acc[i] = BinaryClass(out$`BC Test`)$Accuracy
}

plot(acc, type = 'b', main = 'Accuracy Comparison', xlab = 'Clusters', ylab = 'Acc')
}
}
