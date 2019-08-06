#' Seek the pathway of a gen or genes
#'
#'This function return a data.frame with the pathways associated with the input gen,
#'the data.frame contains the ID, name of pathway and p value
#'
#'
#' @param x A gen symbol in character
#'
#'
#' @return
#'A data.frame with the pathways of the gen from reactome.org
#'
#' @export
#'
#' @examples
#'MAPT <- seeker_gen_pathway("MAPT")
#'

seeker_gen_pathway <- function(x) {
  library(jsonlite)
  server="https://reactome.org/AnalysisService/identifier/"
  pValue_Reactome= list()
  name_Reactome= list()
  pathID_Reactome = list()
  informacion_Reactome <- paste(x, "/projection", sep = "", collapse = NULL)
  url_reactome <- file.path(server,informacion_Reactome, sep = "")
  datos <- fromJSON(url_reactome)
  paths<-datos[["pathways"]]
  paths_select <- data.frame(ID=paths$stId,
                             Path_name=paths$name,
                             pvalue=paths$entities$pValue)
  return(paths_select)
}




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
