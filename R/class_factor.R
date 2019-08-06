#' Class factor
#'
#'seek_gen_pathway for the class factor
#'
#' @param x A genes symbols in class factor
#'
#' @return
#' A data.frame with the pathways of the gen from reactome.org
#' @export
#'
#' @importFrom
#' jsonlite fromJSON
#'
#'
#' @examples
#'seeker_gen_pathway(factor("MAPT"))
#'@export
seeker_gen_pathway.factor <- function(x) {

  message(paste(Sys.time(), 'empezando a correr `seeker_gen_pathway`'))
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
}
