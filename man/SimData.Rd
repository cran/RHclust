\name{SimData}
\alias{SimData}
\title{GE, CPG, SNP Simulated Data}
\description{
Simulated data generator containing continuous variables representing gene expression (GE) data and DNA methylation data as M-values (GPG), and  categorical variable representing single nucleotide polymorphisms (SNP). GE and CPG data are simulated from a normal distribution and SNP data is simulated from a multinomial distribution.
}
\usage{
SimData(seed = NULL, gene = 36,
        k = c(33,33,34),
        GEbar = 5, GEsd = 0.5,
        CPGbar = 4, CPGsd = 0.5,
        SameCPG = FALSE, SameSNP = FALSE,
        ProbDist = NULL)
}
\arguments{
\item{seed}{Set specified seed for reproducibility}
\item{gene}{Numeric input that specifies the number genes}
\item{k}{Cluster pattern/distribution across subjects formatted as a vector, i.e. c(33,33,34) representing 33 subjects in the first cluster, 33 in the second cluster, and 34 in the third cluster.}
\item{GEbar}{Optional numeric input to change the mean distribution of GE data}
\item{GEsd}{Optional numeric input to change the standard deviation of GE data}
\item{CPGbar}{Optional numeric input to change the mean distribution of CPG data}
\item{CPGsd}{Optional numeric input to change the standard deviation of CPG data}
\item{SameCPG}{Logical value that if set to True sets the distribution of each CPG cluster around the same mean}
\item{SameSNP}{Logical value that if set to True changes the probability distribution of SNPs to be the same per cluster}
\item{ProbDist}{Optional list input that allows the change of SNP probability distributions per cluster. Default list stops at 5 cluster distributions. Default ProbDist = list(c(0.50,0.25,0.25), c(0.20,0.55,0.25), c(0.30,0.15,0.55), c(0.20,0.50,0.30), c(0.45,0.20,0.35))}
}

\details{
SimData simply creates simulated data that aims to represent real world data for gene expression (GE), DNA methylations (CPG),
and single neucleotide polymorphisms (SNP). The goal of this function is to allow the user the ability to manipulate their data
for testing of the main VIP() function.
}

\value{
\item{Clusters}{Vector of cluster assignment for each subject.}

\item{Vec}{Numeric representation of values per cluster used for sensitivity measures.}

\item{GE}{Simulated continuous data for GE. Means of each cluster changes by a factor of 5 with default standard deviation of 0.5.}

\item{CPG}{Simulated continuous data for CPG. Means of each cluster changes by a factor of 4 with default standard deviation of 0.5.}

\item{SNP}{Simulated categorical data for SNP.}

\item{GE_Index}{Index names for GE.}

\item{CPG_Index}{Index names for CPG.}

\item{SNP_Index}{Index names for SNP.}

}

\author{
jkhndwrk@memphis.edu
}

\examples{
\donttest{
# Generating simulated data
sd = SimData()

## Specifying seed, genes, and clusters
# sd = SimData(seed = 42, gene = 18, c(10,40,50))

# Specifying probability distribution of SNP data
l = list( c(10, 35, 55),
          c(60, 20, 20),
          c(25, 50, 25))
sd = SimData(1, g = 36, ProbDist = l)

sd = SimData(1, g = 36, c(33,33,34))

# Same SNP distribution across 2 clusters
ssnp = SimData(g = 36, k = c(10,90), SameSNP = TRUE)

# Same CpG distribution across 2 cluster
scpg = SimData(g = 36, k = c(10,90), SameCPG = TRUE)
}
}
