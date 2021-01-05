#' seeker Single Nucleotide Polymorphism consequences
#'
#' @param ID A Single Nucleotide Polymorphism ID ("rs") in character or data.frame
#'
#' @return
#' A list with the consequences based on a variant identifier
#'
#' @source
#' https://rest.ensembl.org
#'
#' @author
#' Erick Cuevas-Fernández
#'
#' @importFrom
#' jsonlite fromJSON
#'
#' @importFrom
#' purrr transpose safely
#'
#' @importFrom
#' furrr future_map
#'
#' @importFrom
#' future plan multiprocess
#'
#'
#' @examples
#' seeker_snp_consequences("rs243034")
#'
#' @rdname seeker_snp_consequences
#' @export seeker_snp_consequences
seeker_snp_consequences <- function(ID) {
  UseMethod("seeker_snp_consequences")
}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_consequences
#' @export
seeker_snp_consequences.character <- function(ID){
  if(length(ID)==1){
    server <- "https://rest.ensembl.org/vep/human/id/"
    ligas <- paste0(server,ID,"?")

    response <- jsonlite::fromJSON(ligas)
    return(response)
  } else{
    df <- data.frame(gene = ID)
    con_result <- seeker_snp_consequences(df)
    return(con_result)
  }

}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_consequences
#' @export
seeker_snp_consequences.factor <- function(ID){
  if(length(ID)==1){
    server <- "https://rest.ensembl.org/vep/human/id/"
    ligas <- paste0(server,ID,"?")

    response <- jsonlite::fromJSON(ligas)
    return(response)
  } else{
    df <- data.frame(gene = ID)
    con_result <- seeker_snp_consequences(df)
    return(con_result)
  }
}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_consequences
#' @export
seeker_snp_consequences.data.frame <- function(ID){

  ID1 <- as.matrix(ID)
  server <- "https://rest.ensembl.org/vep/human/id/"
  ligas <- paste0(server,ID1,"?")
  future::plan("multicore")
  contents <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                .progress = FALSE)
  contents_1 <- purrr::transpose(contents)
  contents_request_first <- contents_1[["result"]]
  # contents_request_first[sapply(contents_request_first, is.null)] <- NULL

  ID1 <- ID1[sapply(contents_request_first, is.null)]
  ligas <- paste0(server,ID1,"?")
  future::plan("multicore")
  contents <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                .progress = FALSE)
  contents_1 <- purrr::transpose(contents)
  contents_request_second <- contents_1[["result"]]

  contents_request <- c(contents_request_first, contents_request_second)
  contents_request[sapply(contents_request, is.null)] <- NULL

  return(contents_request)

}
