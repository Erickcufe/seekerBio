#' seeker Single Nucleotide Polymorphism Linkage-Desequilibrium
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
#' @importFrom
#' jsonlite fromJSON
#'
#' @author
#' Erick Cuevas-Fernández
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

  server_2 <- "https://rest.ensembl.org/ld/human/"
  link1 <- "?window_size="
  link2 <- ";d_prime="
  link3 <- ";content-type=application/json"

  URL_LD <- paste0(server_2, snp,"/",population, link1, window_size,
                   link2, d_prime, link3)

  r <- fromJSON(URL_LD)
}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_ld
#' @export
seeker_snp_ld.data.frame <- function(ID, population = "1000GENOMES:phase_3:MXL",
                                     window_size = 500, d_prime = 0) {

  message(paste(Sys.time(), 'Running `seeker_snp_ld` for data.frame'))

  ID1 <- as.matrix(ID)
  server_2 <- "https://rest.ensembl.org/ld/human/"
  link1 <- "?window_size="
  link2 <- ";d_prime="
  link3 <- ";content-type=application/json"
  URL_LD <- paste0(server_2, ID1,"/",population, link1, window_size,
                   link2, d_prime, link3)

  mydf <- ID[NULL,]
  for(i in 1:length(URL_LD)){

    pop_result <- fromJSON(URL_LD[i])

    mydf <- rbind(mydf, pop_result)
  }
  return(mydf)
}

