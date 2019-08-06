#' Seek the pathway of a gen or genes
#'
#'This function return a data.frame with the pathways associated with the input gen,
#'the data.frame contains the ID, name of pathway and p value
#'
#'
#' @param x A gen symbol in character or more that one gen in a data.frame
#'
#'
#' @return
#'A data.frame with the pathways of the gen from reactome.org
#'
#' @export
#'
#' @import
#' jsonlite
#'
#'@author
#'Cuevas-Fernandez Erick
#'
#'
#' @examples
#'MAPT <- seeker_gen_pathway("MAPT")
#'df <- data.frame(gen=c("MAPT", "APOE", "MMP12"))
#'seeker_gen_pathway(df)

seeker_gen_pathway <- function(x) {
  UseMethod("seeker_gen_pathway")

}



seeker_gen_pathway.character <- function(x) {

  message(paste(Sys.time(), 'empezando a correr `seeker_gen_pathway`'))
  server="https://reactome.org/AnalysisService/identifier/"
  pValue_Reactome= list()
  name_Reactome= list()
  pathID_Reactome = list()
  informacion_Reactome <- paste(x, "/projection", sep = "", collapse = NULL)
  url_reactome <- file.path(server,informacion_Reactome, sep = "")
  datos <- fromJSON(url_reactome)
  paths<-datos[["pathways"]]
  paths_select <- data.frame(Gen = rep(x ,length(paths$stId)),
                             ID=paths$stId,
                             Path_name=paths$name,
                             pvalue=paths$entities$pValue)
  return(paths_select)
}

seeker_gen_pathway.factor <- function(x) {

  message(paste(Sys.time(), 'empezando a correr `seeker_gen_pathway`'))
  server="https://reactome.org/AnalysisService/identifier/"
  pValue_Reactome= list()
  name_Reactome= list()
  pathID_Reactome = list()
  informacion_Reactome <- paste(x, "/projection", sep = "", collapse = NULL)
  url_reactome <- file.path(server,informacion_Reactome, sep = "")
  datos <- fromJSON(url_reactome)
  paths<-datos[["pathways"]]
  paths_select <- data.frame(Gen = rep(x ,length(paths$stId)),
                             ID=paths$stId,
                             Path_name=paths$name,
                             pvalue=paths$entities$pValue)
  return(paths_select)
}



seeker_gen_pathway.data.frame <- function(x) {

  message(paste(Sys.time(), 'empezando a correr `seeker_gen_pathway`'))
  # mydf <- data.frame()
  mydf <- x[NULL,]
  for (i in seq_len(nrow(x))) {

  server="https://reactome.org/AnalysisService/identifier/"
  pValue_Reactome= list()
  name_Reactome= list()
  pathID_Reactome = list()
  informacion_Reactome <- paste(x[i,], "/projection", sep = "", collapse = NULL)
  url_reactome <- file.path(server,informacion_Reactome, sep = "")
  datos <- fromJSON(url_reactome)
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


seeker_gen_pathway.default <- function(x) {
  stop(
    "Don't know how to make seeker_gen_pathway <",
    class(x)[[1]], ">",
    call. = FALSE
  )
}



# df <- data.frame(gen=c("MAPT", "APOE", "MMP12"))
# a<- seeker_gen_pathway(df)
# salio <- scrad_gen_pathway("APOE")
# salio[which.min(salio$pvalue),]
#
# .x<-df[1,]
#
# map_dfr(df, scrad_gen_pathway~(.x))

# to prube
# df<-data.frame(Gene=c("MAPT","APOE","MMP12"))
# df_1 <- as.character(df$Gene)
# Gen<-as.matrix(df_1)
# library(jsonlite)
