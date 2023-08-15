\name{VIPcov}
\alias{VIPcov}
\title{Vector in Partition with covariates}
\description{
Clustering of subjects based on similar patterns of gene expression, DNA methylation, and SNPs weighted by covariates.
}
\usage{
VIPcov(Simulated = NULL, SNP = NULL, CPG = NULL, GE = NULL,
                SNPname = NULL, CPGname = NULL, GEname = NULL,
                covariates = NULL, v, lambda = NULL, gamma = 1,
                optimize = c('off','min','slope','elbow'),
                iter_max = 1000, nstart = 5, fit = c('aic','bic'),
                seed = NULL, type = c('Default','NoCPG','NoSNP'),
                ct = c('mean','median'), verbose = FALSE)
}

\arguments{
\item{Simulated}{set to name of simulated data built from SimData(), else set to NULL for real data.}

\item{SNP}{Data frame or data matrix containing categorical SNP data. Input must be in form of N x M, with N rows of subjects and M columns of SNPs. Rownames are permitted. Run SimData()$SNP for examples.}

\item{CPG}{Data frame or data matrix containing numeric CPG data. Input must be in form of N x M, with N rows of subjects and M columns of CPG. Rownames are permitted. Run SimData()$CPG for examples.}

\item{GE}{Data frame or data matrix containing numeric GE data. Input must be in form of N x M, with N rows of subjects and M columns of GE. Rownames are permitted. Run SimData()$GE for examples.}

\item{SNPname}{Names for SNP data. Data must be a data frame of Nx2 dimensions with SNP sites as column 1, and GE indexes in column 2. Order of SNPs must match the order of the SNP columns in the argument SNP. See SimData()$SNP_Index for examples.}

\item{CPGname}{Names for CPG data. Data must be a data frame of Nx2 dimensions with CPG sites as column 1, and GE indexes in column 2. Order of CPGs must match the order of the CPG columns in the argument GE. See SimData()$CPG_Index for examples.}

\item{GEname}{Names for GE data. Data must be a data frame of Nx2 dimensions with GE sites as column 1, and GE indexes in column 2. Order of GEs must match the order of the GE columns in the argument GE. See SimData()$GE_Index for examples.}

\item{covariates}{Data frame or matrix containing numeric and or categorical variables used for weighting the distance matrix of the genetics data.}

\item{v}{Numeric scalar or vector of number for clusters, or a range of clusters with format c(l,u) for cluster l:u}

\item{lambda}{Weighting factor for the covariate data. The default value is calculated by taking the average variance of the numeric variables divided by the average concentraion of the factors for the categorical variables. When non-mixed type data are used, the value is set to 1 unless specified otherwise.}

\item{gamma}{Weighting factor for the genetics distance measure. Default is set to 1}

\item{optimize}{Returned the optimal number of clusters. Input 'min' returns cluster assignment with lowest WSS for clusters in v. Input 'slope' indicates whether the algorithm should pick the lowest WSS value based on the first increasing slope. Input 'elbow' fits a line between the first and last fitted WSS and finds the corresponding cluster with the maximum distance to that line. All but 'slope' return plots.}

\item{iter_max}{Maximum number of iterations allowed.}

\item{nstart}{If nstart > 1, repetitive computations with random initializations are computed and the result with minimum tot_dist is returned.}

\item{fit}{Penalizing factor for WSS of clusters. Can be set to either 'aic' or 'bic'.}

\item{seed}{Optional input to sample the same initial cluster centers.}

\item{type}{Optional input for special cases for data without CPGs or SNP inputs. Options include "Default", "NoSNP", or "NoCPG"}

\item{ct}{Central tendency option for cluster assignment. Options include 'mean' or 'median'.}

\item{verbose}{Logical whether information about the cluster procedure should be given.}
}

\details{
This method is an extension of the VIP() function and allows for the clustering method to be further weighted by the covarite distance. The covariate distance measure is similar to that seen in \emph{kproto()} where the euclidian measure of the numerics data is weighted by matching categorical variables and a weighting factor. Clsuters are assigned based on these two distance measures where the Central tendancy of numeric data can be set to either mean or median with input \emph{ct}.

Data must be ordered such that rows in each data set correspond to the same subject and order of the indexes match the order of the columns in the data. The current algorithm does not allow for any missing data. The aim is for \emph{GE}, \emph{CPG}, and \emph{SNP} data to be clustered into \emph{v} groups after being weighted by covariates such that within sum of squares is minimized. If groups of clusters are close, the algorithm may not converge correctly and signals a warning if cluster size is reduced.

Optimization functionality was used for simulated data analysis, but is allowed for user exploratory analysis as well. \emph{'min'} simply returns the lowest fitted WSS \emph{fit} parameter. \emph{'slope'} loops through clusters in \emph{v} and returns the cluster based on the first increasing slope of fitted WSS. For example, if AIC output is c(100,80,35,50), cluster 3 would be returned since the slope increases from 3 to 4. If there is no increasing slope, the \emph{'min'} optimizer will be returned. \emph{'elbow'} seeks to find the elbow of the plot based on saturation point. This worked the best for simulation studies but requires more clusters to make proper predictions, in our case it required a range of at least 5 clusters c(1,5) to search to correctly identify the 3 simulated clusters. For ease of exploratory analysis, \emph{v}=1 is allowed. \emph{'lambda'} is used as a weighting factor for the covariate distance measure, while \emph{'gamma'} is the weighting factor for the genetics distance.
}

\value{
\item{size}{Number of subjects assigned to each cluster.}

\item{cluster}{Vector of cluster assignment.}

\item{GECenters}{Matrix of cluster centers for GE.}

\item{CPGCenters}{Matrix of cluster centers for CPG.}

\item{SNPCenters}{Matrix of cluster centers for SNP.}

\item{Adjusted}{Logical if the genetics distance measure is weighted by covariates.}

\item{Lambda}{Lambda value used to weight the covariate distance measure.}

\item{within}{Vector of within cluster sum of squares with one component per cluster.}

\item{tot_within}{Sumed total of within-cluster sum of squares.}

\item{Moved}{Number of iterations before convergence.}

\item{AIC}{Value of tot_within with aic penalizer.}

\item{BIC}{Value of tot_within with bic penalizer.}

\item{outputPlot}{Returns the tot_within, aic, bic, and v values for ploting.}
}

\author{
jkhndwrk@memphis.edu
}

\references{
Hartigan, J. A. and Wong, M. A. (1979). Algorithm AS 136: A K-means clustering algorithm. Applied Statistics, 28, 100--108. 10.2307/2346830.

Z. Huang. (1998) Extensions to the k-means algorithm for clustering large data sets with categorical variables.
Data Mining and Knowledge Discovery, 2:283–304, doi: 10.1023/A:1009769707641.

Szepannek G (2018). “clustMixType: User-Friendly Clustering of Mixed-Type Data in R.” The R Journal, 200-208. doi:10.32614/RJ-2018-048, https://doi.org/10.32614/RJ-2018-048.
}

\examples{
\donttest{
## simple output of 3 clusters assignments
sd = SimData(1, g = 36, c(33,33,34))
VIPout = VIPcov(sd, v = 3, covariates = sd$Covariates)

# loop through clusters 1-10 and outputs plot of WSS, AIC, and BIC
VIPout = VIPcov(sd, v = c(1,10), covariates = sd$Covariates)

# loop through clusters 1-10 but picks first instance of increasing slope
VIPout = VIPcov(sd, v = c(1,10), optimize = 'slope', covariates = sd$Covariates)


# Individual inputs
sd = SimData(1, g = 36, k = c(33,33,34))
VIPout = VIPcov(SNP = sd$SNP, CPG = sd$CPG, GE = sd$GE,
            SNPname = sd$SNP_Index, CPGname = sd$CPG_Index,
            GEname = sd$GE_Index,
            v = c(1,5), optimize = 'off', nstart = 5, covariates = sd$Covariates)


## Varying clusters
sd = SimData(k = c(10,40,50))
out = VIPcov(sd, v = c(1,6), optimize = 'elbow', nstart = 30, covariates = sd$Covariates)

}
}
