#' Class data.frame
#'
#'seek_gen_pathway for the class data.frame
#'
#' @param x A genes symbols in dataframe
#'
#' @return
#' A data.frame with the pathways of the gen from reactome.org
#' @export
#'
#'
#'
#' @importFrom
#' jsonlite fromJSON
#'
#'
#' @examples
#' df <- data.frame(gen=c("MAPT", "APOE", "MMP12"))
#'seeker_gen_pathway(df)
seeker_gen_pathway.data.frame <- function(x) {

  message(paste(Sys.time(), 'empezando a correr `seeker_gen_pathway`'))
  mydf <- x[NULL,]
  for (i in seq_len(nrow(x))) {

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


