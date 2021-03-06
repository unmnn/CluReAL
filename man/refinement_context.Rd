% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refinement-context.R
\name{refinement_context}
\alias{refinement_context}
\title{Create the refinement context}
\usage{
refinement_context(x, y, cc, gv)
}
\arguments{
\item{x}{The input data (numeric matrix).}

\item{y}{Cluster membership indices (integer vector). -1 indicates noise.}

\item{cc}{The \code{\link[cluster_context]{cluster context}}.}

\item{gv}{The \code{\link[gval]{GOI cluster validity}} indices object.}
}
\value{
A list with four elements:
\itemize{
\item \code{mm}: Cluster multimodality flag (logical vector of length k).
\item \code{k_dens}: Cluster-individual relative density (double vector of
length k).
\item \code{global_c_dens}: Global density (double).
\item \code{kinship}: cluster kinship matrix (k x k matrix): 0-itself,
1-parent and child, 2-relatives, 3-close friends, 4-acquaintances,
5-unrelated.
}
}
\description{
Create the refinement context
}
\examples{
# TODO
}
