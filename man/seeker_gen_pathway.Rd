% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seeker_gen_pathway.R
\name{seeker_gen_pathway}
\alias{seeker_gen_pathway}
\alias{seeker_gen_pathway.character}
\alias{seeker_gen_pathway.factor}
\alias{seeker_gen_pathway.data.frame}
\alias{seeker_gen_pathway.default}
\title{Seeker the pathway of a gen or genes}
\source{
https://reactome.org
}
\usage{
seeker_gen_pathway(x)

\method{seeker_gen_pathway}{character}(x)

\method{seeker_gen_pathway}{factor}(x)

\method{seeker_gen_pathway}{data.frame}(x)

\method{seeker_gen_pathway}{default}(x)
}
\arguments{
\item{x}{A gen symbol}
}
\value{
A data.frame with the pathways of the gen from reactome.org

\code{NULL}

\code{NULL}

\code{NULL}

\code{NULL}
}
\description{
seeker_gene_pathway is a generic function to produce a data.frame with the pathways associated with the input gen,
the data.frame contains the ID, name of pathway and p value. The function invokes particular methods wich depend on the class of the first argument
}
\examples{
MAPT <- seeker_gen_pathway("MAPT")
df <- data.frame(gen=c("MAPT", "APOE", "MMP12"))
seeker_gen_pathway(df)
}
\author{
Erick Cuevas Fernández

Heriberto Manuel Rivera
}
