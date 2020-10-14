#' seeker Single Nucleotide Polymorphism frequency
#'
#' A generic function to search the population frequency of a specific study
#'
#' @param ID A Single Nucleotide Polymorphism ID ("rs") in character or data.frame
#' @param study A study of population frequency. Default ("1000GENOMES:phase3")
#'
#' @return
#' A data.frame with the allel frequency in all the population of the given study
#'
#' @source
#' https://rest.ensembl.org
#'
#' @author
#' Erick Cuevas-Fernández
#'
#' Heriberto Manuel Rivera
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
#' @examples
#' seeker_snp_freq("rs56116432")
#'
#' df <- data.frame(c("rs56116432","rs10878307", "rs7133914", "rs11564148", "rs3761863", "rs10878245"))
#' seeker_snp_freq(df)
#'
#' @rdname seeker_snp_freq
#' @export seeker_snp_freq
seeker_snp_freq <- function(ID, study = "1000GENOMES:phase_3") {
  UseMethod("seeker_snp_freq")
}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_freq
#' @export
seeker_snp_freq.character <- function(ID, study = "1000GENOMES:phase_3"){
  # message(paste(Sys.time(), 'Running `seeker_snp_freq` for character'))

  if (length(ID)==1){
  server <- "http://rest.ensembl.org/variation/human/"
  ligas <- paste0(server, ID,"?pops=1;content-type=application/json")

    r <- fromJSON(ligas)
    pop <- r[["populations"]]

    seleccion <- stringr::str_detect(pop$population, study)
    SNP <- c(rep(ID, length(pop[seleccion,])))
    pop_result <- cbind(SNP = ID, pop[seleccion,])
    pop_result$submission_id <- NULL

  return(pop_result)
  } else {

    df <- data.frame(gene = ID)
    pop_result <- seeker_snp_freq(df)
    return(pop_result)

  }
}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_freq
#' @export
seeker_snp_freq.factor <- function(ID, study = "1000GENOMES:phase_3"){
  # message(paste(Sys.time(), 'Running `seeker_snp_freq` for factor'))

  if (length(ID)==1){
    server <- "http://rest.ensembl.org/variation/human/"
    ligas <- paste0(server, ID,"?pops=1;content-type=application/json")

    r <- fromJSON(ligas)
    pop <- r[["populations"]]

    seleccion <- stringr::str_detect(pop$population, study)
    SNP <- c(rep(ID, length(pop[seleccion,])))
    pop_result <- cbind(SNP = ID, pop[seleccion,])
    pop_result$submission_id <- NULL

    return(pop_result)
  } else {

    df <- data.frame(gene = ID)
    pop_result <- seeker_snp_freq(df)
    return(pop_result)

  }
}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_freq
#' @export
seeker_snp_freq.data.frame <- function(ID, study = "1000GENOMES:phase_3"){
  # message(paste(Sys.time(), 'Running `seeker_snp_freq` for data.frame'))
  ID <- unique(ID)
  ID1 <- as.matrix(ID)
  server <- "http://rest.ensembl.org/variation/human/"
  ligas <- paste0(server, ID1,"?pops=1;content-type=application/json")
  future::plan("multiprocess")
  contents <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                .progress = FALSE)
  contents_1 <- purrr::transpose(contents)
  while(sum(!sapply(contents_1[["error"]], is.null)) == length(contents_1[["error"]])){
    future::plan("multiprocess")
    contents <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                  .progress = FALSE)
    contents_1 <- purrr::transpose(contents)
    # message(contents_1[["error"]][[1]])
  }
  contents_request_first <- contents_1[["result"]]

  if(sum(!sapply(contents_1[["error"]], is.null)) == 0){
    contents_request <- contents_1[["result"]]
    mydf <- data.frame()
    for(i in 1:length(contents_request)){
      contents_request <- contents_request_first
      pop <- contents_request[[i]][["populations"]]
      if (length(pop)==0){
        next()
      }
      seleccion <- stringr::str_detect(pop$population, study)
      if (!is.null(pop) & length(pop[seleccion,1])!=0){
        pop_select <- pop[seleccion,]
        pop_result <- data.frame(SNP = contents_request[[i]]$name,
                                 pop_select)
        pop_result$submission_id <- NULL
        mydf <- rbind(mydf, pop_result)
      } else{
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
      future::plan("multiprocess")
      contents <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                    .progress = FALSE)
      contents_1 <- purrr::transpose(contents)
      # message(contents_1[["error"]][[1]])
    }
    contents_request_second <- contents_1[["result"]]
    ID3 <- ID2[sapply(contents_request_second, is.null)]
    if(length(ID3) > 1){
      ligas <- paste0(server, ID3,"?pops=1;content-type=application/json")
      future::plan("multiprocess")
      contents_2 <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                      .progress = FALSE)
      contents_3 <- purrr::transpose(contents_2)
      while(sum(!sapply(contents_3[["error"]], is.null)) > 0){
        future::plan("multiprocess")
        contents_2 <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                      .progress = FALSE)
        contents_3 <- purrr::transpose(contents_2)
        # message(contents_3[["error"]][[1]])
      }
      contents_3_request <-  contents_3[["result"]]
      contents_request <- c(contents_request_first, contents_request_second,
                            contents_3_request)
    } else{
      contents_request <- c(contents_request_first, contents_request_second)
    }

    contents_request[sapply(contents_request, is.null)] <- NULL
    mydf <- data.frame()
    for(i in 1:length(contents_request)){
      pop <- contents_request[[i]][["populations"]]
      if (length(pop)==0){
        next()
      }
      seleccion <- stringr::str_detect(pop$population, study)
      if (!is.null(pop) & length(pop[seleccion,1])!=0){
        pop_select <- pop[seleccion,]
        pop_result <- data.frame(SNP = contents_request[[i]]$name,
                                 pop_select)
        pop_result$submission_id <- NULL
        mydf <- rbind(mydf, pop_result)
      } else{
        next()
      }
    }
  }
  return(mydf)
}



