#' seeker Single Nucleotide Polymorphism (SNP) arquitecture
#'
#' seeker_snp_arq is a generic function to obtain the genome arquitecture data of a given SNP
#'
#' @param ID A Single Nucleotide Polymorphism (SNP) ID "rs" in character or data.frame
#'
#' @return
#' A data.frame with the genomic arquitecture
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
#' @author
#' Erick Cuevas-Fernández
#'
#' Heriberto Manuel Rivera
#'
#' @source
#' https://rest.ensembl.org
#'
#' @examples
#' seeker_snp_arq("rs56116432")
#'
#' df <- data.frame(c("rs56116432","rs10878307", "rs7133914", "rs11564148", "rs3761863", "rs10878245"))
#' seeker_snp_arq(df)
#'
#' @rdname seeker_snp_arq
#' @export seeker_snp_arq
seeker_snp_arq <- function(ID){
  UseMethod("seeker_snp_arq")
}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_arch
#' @export
seeker_snp_arq.character <- function(ID){

  # message(paste(Sys.time(), 'Running `seeker_snp_arch` for character'))

  if (length(ID)==1){
    server <- "http://rest.ensembl.org/variation/human/"
    ligas <- paste0(server, ID,"?pops=1;content-type=application/json")

    r <- fromJSON(ligas)
    pop <- r[["mappings"]]
    pop_result <- data.frame(SNP = ID, pop)

    return(pop_result)
  } else {

    df <- data.frame(gene = ID)
    pop_result <- seeker_snp_arq(df)
    return(pop_result)

  }

}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_arch
#' @export
seeker_snp_arq.factor <- function(ID){

  # message(paste(Sys.time(), 'Running `seeker_snp_arch` for factor'))

  if (length(ID)==1){
    server <- "http://rest.ensembl.org/variation/human/"
    ligas <- paste0(server, ID,"?pops=1;content-type=application/json")

    r <- fromJSON(ligas)
    pop <- r[["mappings"]]
    pop_result <- data.frame(SNP = ID, pop)

    return(pop_result)
  } else {

    df <- data.frame(gene = ID)
    pop_result <- seeker_snp_arq(df)
    return(pop_result)

  }

}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_arch
#' @export
seeker_snp_arq.data.frame <- function(ID){
  # message(paste(Sys.time(), 'Running `seeker_snp_arch` for data.frame'))
  ID <- unique(ID)
  ID1 <- as.matrix(ID)
  server <- "http://rest.ensembl.org/variation/human/"
  ligas <- paste0(server, ID1,"?pops=1;content-type=application/json")
  future::plan("multiprocess")
  contents <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                .progress = FALSE)
  contents_1 <- purrr::transpose(contents)
  while(sum(!sapply(contents_1[["error"]], is.null)) == length(contents_1[["error"]])){
    contents <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                  .progress = FALSE)
    contents_1 <- purrr::transpose(contents)
    message(contents_1[["error"]][[1]])
  }
  contents_request_first <- contents_1[["result"]]

  if(sum(!sapply(contents_1[["error"]], is.null)) == 0){
    contents_request <- contents_1[["result"]]
    mydf <- data.frame()
    for (i in 1:length(contents_request)){
      pop <- contents_request[[i]][["mappings"]]
      if (length(pop)==0){
        next()
      }
      if(!is.null(pop) & length(pop[,1])!=0){
        pop_result <- data.frame(SNP = contents_request[[i]]$name, pop)
        mydf <- rbind(mydf, pop_result)
      } else {
        next()
      }
    }
  } else {
    ID2 <- ID1[sapply(contents_request_first, is.null)]
    server <- "http://rest.ensembl.org/variation/human/"
    ligas <- paste0(server, ID2,"?pops=1;content-type=application/json")
    future::plan(multiprocess)
    contents <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                  .progress = FALSE)
    contents_1 <- purrr::transpose(contents)
    while(sum(!sapply(contents_1[["error"]], is.null)) == length(contents_1[["error"]])){
      contents <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                    .progress = FALSE)
      contents_1 <- purrr::transpose(contents)
      message(contents_1[["error"]][[1]])
    }
    while(sum(!sapply(contents_1[["result"]], is.null)) < length(ID2)){
      contents <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                    .progress = FALSE)
      contents_1 <- purrr::transpose(contents)
      message(contents_1[["error"]][[1]])
    }
    contents_request_second <- contents_1[["result"]]
    contents_request <- c(contents_request_first, contents_request_second)
    mydf <- data.frame()
    for (i in 1:length(contents_request)){
      pop <- contents_request[[i]][["mappings"]]
      if (length(pop)==0){
        next()
      }
      if(!is.null(pop) & length(pop[,1])!=0){
        pop_result <- data.frame(SNP = contents_request[[i]]$name, pop)
        mydf <- rbind(mydf, pop_result)
      } else {
        next()
      }
    }
  }

  return(mydf)
}