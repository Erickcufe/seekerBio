#' seeker Single Nucleotide Polymorphism Linkage-Desequilibrium
#'
#' A function that returns LD values between the SNP and all other SNP in a window centered around the given SNP.
#'
#' @param ID A Single Nucleotide Polymorphism ID ("rs") in character or data.frame
#' @param population Population for which to compute LD. Default "1000GENOMES:phase_3:MXL"
#' @param window_size Window size in kb. The maximum allowed value for the window size is 500 kb.
#' LD is computed for the given variant and all variants that are located within the specified window.
#' Default 500
#' @param d_prime Measure of LD. If D' is provided only return pairs of variants whose D' value is equal to or greater than the value provided.
#' Default 0
#'
#' @return
#' A data.frame with the LD information of SNP in a specific population
#'
#' @source
#' https://rest.ensembl.org
#'
#' @importFrom
#' jsonlite fromJSON
#'
#' @importFrom
#' purrr map
#'
#' @author
#' Erick Cuevas-Fernández
#' Heriberto Manuel Rivera
#'
#' @examples
#'seeker_snp_ld("rs56116432")
#'
#'df <- data.frame(c("rs56116432","rs10878307", "rs7133914", "rs11564148", "rs3761863", "rs10878245"))
#'seeker_snp_ld(df)
#'
#' @rdname seeker_snp_ld
#' @export seeker_snp_ld
seeker_snp_ld <- function(ID, population = "1000GENOMES:phase_3:MXL",
                          window_size = 500, d_prime = 0){
  UseMethod("seeker_snp_ld")
}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_ld
#' @export
seeker_snp_ld.character <- function(ID, population = "1000GENOMES:phase_3:MXL",
                          window_size = 500, d_prime = 0){
  message(paste(Sys.time(), 'Running `seeker_snp_ld` for character'))

  if (length(ID)==1){
  server_2 <- "https://rest.ensembl.org/ld/human/"
  link1 <- "?window_size="
  link2 <- ";d_prime="
  link3 <- ";content-type=application/json"

  URL_LD <- paste0(server_2, snp,"/",population, link1, window_size,
                   link2, d_prime, link3)

  r <- fromJSON(URL_LD)
  return(r)
  } else {
    df <- data.frame(ID)
    ld_result <- seeker_snp_ld(df)
    return(ld_result)
  }
}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_ld
#' @export
seeker_snp_ld.factor <- function(ID, population = "1000GENOMES:phase_3:MXL",
                                    window_size = 500, d_prime = 0){
  message(paste(Sys.time(), 'Running `seeker_snp_ld` for factor'))

  if (length(ID)==1){
  server_2 <- "https://rest.ensembl.org/ld/human/"
  link1 <- "?window_size="
  link2 <- ";d_prime="
  link3 <- ";content-type=application/json"

  URL_LD <- paste0(server_2, snp,"/",population, link1, window_size,
                   link2, d_prime, link3)

  r <- fromJSON(URL_LD)
  return(r)
  }  else {
    df <- data.frame(ID)
    ld_result <- seeker_snp_ld(df)
    return(ld_result)
  }
}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_ld
#' @export
seeker_snp_ld.data.frame <- function(ID, population = "1000GENOMES:phase_3:MXL",
                                     window_size = 500, d_prime = 0){

  message(paste(Sys.time(), 'Running `seeker_snp_ld` for data.frame'))

  ID1 <- as.matrix(ID)
  server_2 <- "https://rest.ensembl.org/ld/human/"
  link1 <- "?window_size="
  link2 <- ";d_prime="
  link3 <- ";content-type=application/json"
  URL_LD <- paste0(server_2, ID1,"/",population, link1, window_size,
                   link2, d_prime, link3)

  contents <- purrr::map(URL_LD, safely(jsonlite::fromJSON))
  contents_1 <- purrr::transpose(contents)
  contents_request <- contents_1[["result"]]

  mydf <- data.frame()
  for(i in 1:length(URL_LD)){
    if (!is.null(contents_request[[i]])){
      mydf <- rbind(mydf, contents_request[[i]])
    }
    return(mydf)
  }
}

