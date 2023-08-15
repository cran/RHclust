\name{VIPnoSNP}
\alias{VIPnoSNP}
\title{Vector in Partition without SNP data}
\description{
Clustering of subjects based on similar patterns of gene expression and DNA methylation.
}
\usage{
VIPnoSNP(Simulated = NULL, CPG = NULL, GE = NULL,
                CPGname = NULL, GEname = NULL, v,
                optimize = c('off','min','slope','elbow'),
                iter_max = 1000, nstart = 5, fit = c('aic','bic'),
                seed = NULL, ct = c('mean','median'), verbose = FALSE)
}

\arguments{
\item{Simulated}{set to name of simulated data built from SimData(), else set to NULL for real data.}

\item{CPG}{Data frame or data matrix containing numeric CPG data. Input must be in form of N x M, with N rows of subjects and M columns of CPG. Rownames are permitted. Run SimData()$CPG for examples.}

\item{GE}{Data frame or data matrix containing numeric GE data. Input must be in form of N x M, with N rows of subjects and M columns of GE. Rownames are permitted. Run SimData()$GE for examples.}

\item{CPGname}{Names for CPG data. Data must be a data frame of Nx2 dimensions with CPG sites as column 1, and GE indexes in column 2. Order of CPGs must match the order of the CPG columns in the argument GE. See SimData()$CPG_Index for examples.}

\item{GEname}{Names for GE data. Data must be a data frame of Nx2 dimensions with GE sites as column 1, and GE indexes in column 2. Order of GEs must match the order of the GE columns in the argument GE. See SimData()$GE_Index for examples.}

\item{v}{Numeric scalar or vector of number for clusters, or a range of clusters with format c(l,u) for cluster l:u}

\item{optimize}{Returned the optimal number of clusters. Input 'min' returns cluster assignment with lowest WSS for clusters in v. Input 'slope' indicates whether the algorithm should pick the lowest WSS value based on the first increasing slope. Input 'elbow' fits a line between the first and last fitted WSS and finds the corresponding cluster with the maximum distance to that line. All but 'slope' return plots.}

\item{iter_max}{Maximum number of iterations allowed.}

\item{nstart}{If nstart > 1, repetitive computations with random initializations are computed and the result with minimum tot_dist is returned.}

\item{fit}{Penalizing factor for WSS of clusters. Can be set to either 'aic' or 'bic'.}

\item{seed}{Optional input to sample the same initial cluster centers.}

\item{ct}{Central tendency option for cluster assignment. Options include 'mean' or 'median'.}

\item{verbose}{Logical whether information about the cluster procedure should be given.}
}

\details{
The details are outlined in the main VIP() function. The only difference in this function is the absence of SNP data.
}

\value{
\item{size}{Number of subjects assigned to each cluster.}

\item{cluster}{Vector of cluster assignment.}

\item{GECenters}{Matrix of cluster centers for GE.}

\item{CPGCenters}{Matrix of cluster centers for CPG.}

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
}

\examples{
\donttest{
# No SNP data
sd = SimData()
noSNPout = VIP(sd, v = c(1,5), optimize = 'off', nstart = 30, type = 'NoSNP')

noSNPout = VIPnoSNP(sd, v = c(1,5), optimize = 'off', nstart = 30)
}
}
