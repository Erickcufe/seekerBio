% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seeker_snp_consequences.R
\name{seeker_snp_consequences}
\alias{seeker_snp_consequences}
\alias{seeker_snp_consequences.character}
\alias{seeker_snp_consequences.factor}
\alias{seeker_snp_consequences.data.frame}
\title{seeker Single Nucleotide Polymorphism consequences}
\source{
https://rest.ensembl.org
}
\usage{
seeker_snp_consequences(ID)

\method{seeker_snp_consequences}{character}(ID)

\method{seeker_snp_consequences}{factor}(ID)

\method{seeker_snp_consequences}{data.frame}(ID)
}
\arguments{
\item{ID}{A Single Nucleotide Polymorphism ID ("rs") in character or data.frame}
}
\value{
A list with the consequences based on a variant identifier

\code{NULL}

\code{NULL}

\code{NULL}
}
\description{
seeker Single Nucleotide Polymorphism consequences
}
\examples{
seeker_snp_consequences("rs243034")

}
\author{
Erick Cuevas-Fernández
}
