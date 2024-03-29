% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deltaccd.R
\name{calcDeltaCCD}
\alias{calcDeltaCCD}
\title{Calculate delta clock correlation distance.}
\usage{
calcDeltaCCD(
  refCor,
  emat,
  groupVec,
  groupNormal,
  refEmat = NULL,
  nPerm = 1000,
  geneNames = NULL,
  dopar = FALSE,
  scale = FALSE
)
}
\arguments{
\item{refCor}{Correlation matrix to be used as the reference, such as comes
from \code{\link[=getRefCor]{getRefCor()}}. Should contain Spearman correlation values.}

\item{emat}{Matrix of expression values, where each row corresponds to a gene
and each column corresponds to a sample. The rownames and colnames of
\code{refCor} should be present in the rownames of \code{emat}. For the p-value
calculation, it is important that \code{emat} include all measured genes, not
just those in \code{refCor}.}

\item{groupVec}{Vector indicating the group to which group each sample
belongs. It's ok for groupVec to have more than two groups.}

\item{groupNormal}{Value indicating the group in groupVec that corresponds to
normal or healthy. Other groups will be compared to this group.}

\item{refEmat}{Optional expression matrix for calculating co-expression for
the reference, with the same organization as \code{emat}. Only used if \code{refCor}
is not provided.}

\item{nPerm}{Number of permutations for assessing statistical significance.}

\item{geneNames}{Optional vector indicating a subset of genes in \code{refCor},
\code{emat}, and/or \code{refEmat} to use for calculating the CCD.}

\item{dopar}{Logical indicating whether to process features in parallel. Make
sure to register a parallel backend first.}

\item{scale}{Logical indicating whether to use scaled CCDs to calculate
difference.}
}
\value{
A data.table with columns for group 1, group 2, deltaCCD, and
p-value. In each row, the deltaCCD is the CCD of group 2 minus the CCD of
group 1, so group 1 corresponds to \code{groupNormal}.
}
\description{
Calculate the difference between the clock correlation distances (CCDs),
relative to a reference, for two groups of samples. Statistical significance
is calculated using permutation of the samples that belong to either of those
two groups.
}
\examples{
set.seed(35813)

refCor = getRefCor()
deltaCcdResult = calcDeltaCCD(
  refCor, GSE19188$emat, GSE19188$groupVec, 'healthy', nPerm = 100)

}
\seealso{
\code{\link[=getRefCor]{getRefCor()}}, \code{\link[=calcCCD]{calcCCD()}}, \code{\link[=plotHeatmap]{plotHeatmap()}}
}
