% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seeker_gwas.R
\name{seeker_gwas}
\alias{seeker_gwas}
\alias{seeker_gwas.character}
\title{seeker gwas}
\source{
https://www.ebi.ac.uk/gwas/
}
\usage{
seeker_gwas(trait)

\method{seeker_gwas}{character}(trait)
}
\arguments{
\item{trait}{A keyword to search in GWAS}
}
\value{
A data.frame with all the information from https://www.ebi.ac.uk/gwas/home

\code{NULL}
}
\description{
seeker_gwas is a generic function that allow to download a data.frame of all Genome-wide
Association Study of the input trait
}
\examples{
seeker_gwas("multiple sclerosis")
seeker_gwas("fontotemporal")
}
\author{
Erick Cuevas Fernández

Heriberto Manuel Rivera
}
