#' Seeker the pathway of a gen or genes
#'
#'seeker_gene_pathway is a generic function to produce a data.frame with the pathways associated with the input gen,
#'the data.frame contains the ID, name of pathway and p value. The function invokes particular methods wich depend on the class of the first argument
#'
#'
#' @param x A gen symbol
#'
#'
#' @return
#'A data.frame with the pathways of the gen from reactome.org
#'
#' @importFrom
#'jsonlite fromJSON
#'
#'
#' @author
#'Erick Cuevas-Fernández
#'
#' @source
#' https://reactome.org
#'
#' @examples
#'MAPT <- seeker_gen_pathway("MAPT")
#'df <- data.frame(gen=c("MAPT", "APOE", "MMP12"))
#'seeker_gen_pathway(df)
#' @rdname seeker_gen_pathway
#' @export seeker_gen_pathway
seeker_gen_pathway <- function(x) {
  UseMethod("seeker_gen_pathway")

}

#' @return \code{NULL}
#'
#' @rdname seeker_gen_pathway
#' @export
seeker_gen_pathway.character <- function(x) {

  message(paste(Sys.time(), 'Running `seeker_gen_pathway` for character'))


  if (length(x)==1){

  server="https://reactome.org/AnalysisService/identifier/"
  pValue_Reactome= list()
  name_Reactome= list()
  pathID_Reactome = list()
  informacion_Reactome <- paste(x, "/projection", sep = "", collapse = NULL)
  url_reactome <- file.path(server,informacion_Reactome, sep = "")
  datos <- jsonlite::fromJSON(url_reactome)
  paths<-datos[["pathways"]]
  paths_select <- data.frame(Gen = rep(x ,length(paths$stId)),
                             ID=paths$stId,
                              Path_name=paths$name,
                             pvalue=paths$entities$pValue)
  return(paths_select)
  } else {

    mydf <- data.frame(gene=x)
    mypaths <- seeker_gen_pathway(mydf)
    return(mypaths)
  }

}

#' @return \code{NULL}
#'
#' @rdname seeker_gen_pathway
#' @export
seeker_gen_pathway.factor <- function(x) {
  message(paste(Sys.time(), 'Running `seeker_gen_pathway` for factor'))
  if (length(x)==1){

  server="https://reactome.org/AnalysisService/identifier/"
  pValue_Reactome= list()
  name_Reactome= list()
  pathID_Reactome = list()
  informacion_Reactome <- paste(x, "/projection", sep = "", collapse = NULL)
  url_reactome <- file.path(server,informacion_Reactome, sep = "")
  datos <- jsonlite::fromJSON(url_reactome)
  paths<-datos[["pathways"]]
  paths_select <- data.frame(Gen = rep(x ,length(paths$stId)),
                             ID=paths$stId,
                             Path_name=paths$name,
                             pvalue=paths$entities$pValue)
  return(paths_select)
  } else {

    mydf <- data.frame(gene=x)
    mypaths <- seeker_gen_pathway(mydf)
    return(mypaths)
  }
}


#' @return \code{NULL}
#'
#' @rdname seeker_gen_pathway
#' @export
seeker_gen_pathway.data.frame <- function(x) {

  message(paste(Sys.time(), 'Running `seeker_gen_pathway` for data.frame'))

  mydf <- x[NULL,]
  for (i in seq_len(nrow(x))) {
  if (x[i,]=="" | x[i,]=="NR") {
    next()
  }
  server="https://reactome.org/AnalysisService/identifier/"
  pValue_Reactome= list()
  name_Reactome= list()
  pathID_Reactome = list()
  informacion_Reactome <- paste(x[i,], "/projection", sep = "", collapse = NULL)
  url_reactome <- file.path(server,informacion_Reactome, sep = "")
  datos <- jsonlite::fromJSON(url_reactome)
  paths<-datos[["pathways"]]
  paths_select <- data.frame(Gen = rep(x[i,] ,length(paths$stId)),
                             ID=paths$stId,
                             Path_name=paths$name,
                             pvalue=paths$entities$pValue)
  mydf <- rbind(mydf, paths_select)

  }

  mypaths <- data.frame(mydf)
  colnames(mypaths) <- c("Gen", "ID", "Path_name", "pvalue")
  return(mypaths)
}

#' @return \code{NULL}
#'
#' @rdname seeker_gen_pathway
#' @export
seeker_gen_pathway.default <- function(x) {
  stop(
    "Don't know how to make seeker_gen_pathway <",
    class(x)[[1]], ">",
    call. = FALSE
  )
}



