#' seeker Single Nucleotide Polymorphism Context
#'
#' A generic function to search the type of variation & gene
#'
#' @param SNP A Single Nucleotide Polymorphism ID ("rs") in character or data.frame
#'
#' @return
#' A data.frame with Type of variation
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
#' seeker_snp_context("rs56116432")
#'
#' df <- data.frame(c("rs56116432","rs10878307", "rs7133914", "rs11564148", "rs3761863", "rs10878245"))
#' seeker_snp_context(df)
#'
#' @rdname seeker_snp_context
#' @export seeker_snp_context
seeker_snp_context <- function(SNP) {
  UseMethod("seeker_snp_context")
}

#' @rdname seeker_snp_context
#' @export
seeker_snp_context.factor <- function(SNP){
  # message(paste(Sys.time(), 'Running `seeker_snp_context` for factor'))

  if (length(SNP)==1){

    SNP <- as.character(SNP)
    seeker_snp_context(SNP=SNP)

  } else {

    df <- data.frame(SNP = SNP)
    pop_result <- seeker_snp_context(SNP=df)
    return(pop_result)


  }

}

#' @rdname seeker_snp_context
#' @export
seeker_snp_context.character <- function(SNP){
  # message(paste(Sys.time(), 'Running `seeker_snp_context` for character'))

  if (length(SNP)==1){

    URL_dbSNP <- "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/"
    server <- "http://rest.ensembl.org/variation/human/"
    ligas_context <- paste0(server, SNP,"?pops=1;content-type=application/json")
    ligas <- paste0(URL_dbSNP,SNP)


    # For context
    future::plan(multiprocess)
    contents <- furrr::future_map(ligas_context, purrr::safely(jsonlite::fromJSON),
                                  .progress = TRUE)
    contents_1 <- purrr::transpose(contents)
    contents_request <- contents_1[["result"]]
    context <- contents_request[[1]]$most_severe_consequence
    snps_context <- data.frame(SNP = SNP, CONTEXT = context)

    return(snps_context)

  } else {
    df <- data.frame(SNP = SNP)
    pop_result <- seeker_snp_context(SNP=df)
    return(pop_result)
  }

}


#' @rdname seeker_snp_context
#' @export
seeker_snp_context.data.frame <- function(SNP){

  # message(paste(Sys.time(), 'Running `seeker_snp_context` for data.frame'))
  SNPs <- unique(SNP)
  SNPs <- as.matrix(SNPs)
  URL_dbSNP <- "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/"
  ligas <- paste0(URL_dbSNP,SNPs)
  server <- "http://rest.ensembl.org/variation/human/"
  ligas_context <- paste0(server, SNPs,"?pops=1;content-type=application/json")

  # For context
  future::plan("multiprocess")
  contents <- furrr::future_map(ligas_context, purrr::safely(jsonlite::fromJSON),
                                .progress = FALSE)
  contents_1 <- purrr::transpose(contents)
  while(sum(!sapply(contents_1[["error"]], is.null)) == length(contents_1[["error"]])){
    future::plan("multiprocess")
    contents <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                  .progress = FALSE)
    contents_1 <- purrr::transpose(contents)
    # message(contents_1[["error"]][[1]])
  }
  contents_request <- contents_1[["result"]]
  ID3 <- SNPs[sapply(contents_request, is.null)]
    if(length(ID3) > 1){
      ligas <- paste0(server, ID3,"?pops=1;content-type=application/json")
      future::plan("multiprocess")
      contents_2 <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                      .progress = FALSE)
      contents_3 <- purrr::transpose(contents_2)
      if(sum(!sapply(contents_3[["error"]], is.null)) == length(contents_3[["error"]])){

        while(sum(!sapply(contents_3[["error"]], is.null)) == length(contents_3[["error"]])){
          contents_2 <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                          .progress = FALSE)
          contents_3 <- purrr::transpose(contents_2)
          error_400 <- vector()
          contents_3[sapply(contents_3[["error"]], is.null)] <- NULL
          for(i in 1:length(contents_1[["error"]])){
            error_400 <- c(error_400,contents_1[["error"]][[i]][["message"]] == "HTTP error 400.")
          }
          if(sum(error_400) >= 1){
            break
          }
        }
      }
    contents_3_request <-  contents_3[["result"]]
    contents_request <- c(contents_request,
                          contents_3_request)
  } else{
    contents_request <- contents_request
  }
  snps_context <- data.frame()
  for (i in 1:length(contents_request)) {
    context <- contents_request[[i]]$most_severe_consequence
    snp <- contents_request[[i]]$name
    if(is.null(context)){
      context <- "none"
      next
    }
    df <- data.frame(SNP = snp, CONTEXT = context)
    snps_context <- rbind(snps_context, df)

  }

  return(snps_context)
}

