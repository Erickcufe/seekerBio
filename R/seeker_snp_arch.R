#' seeker Single Nucleotide Polymorphism (SNP) arquitecture
#'
#' seeker_snp_arch is a generic function to obtain the genome arquitecture data of a given SNP
#'
#' @param ID A Single Nucleotide Polymorphism (SNP) ID "rs" in character or data.frame
#'
#' @return
#' A data.frame with the genomic arquitecture
#'
#' @importFrom
#' jsonlite fromJSON
#'
#' @source
#' https://rest.ensembl.org
#'
#' @examples
#' seeker_snp_arch("rs56116432")
#'
#' df <- data.frame(c("rs56116432","rs10878307", "rs7133914", "rs11564148", "rs3761863", "rs10878245"))
#' seeker_snp_arch(df)
#'
#' @rdname seeker_snp_arch
#' @export seeker_snp_arch
seeker_snp_arch <- function(ID){
  UseMethod("seeker_snp_arch")
}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_arch
#' @export
seeker_snp_arch.character <- function(ID){

  message(paste(Sys.time(), 'Running `seeker_snp_arch` for character'))

  if (length(ID)==1){
    server <- "http://rest.ensembl.org/variation/human/"
    ligas <- paste0(server, ID,"?pops=1;content-type=application/json")

    r <- fromJSON(ligas)
    pop <- r[["mappings"]]
    pop_result <- data_frame(SNP = ID, pop)

    return(pop_result)
  } else {

    df <- data.frame(gene = ID)
    pop_result <- seeker_snp_arch(df)
    return(pop_result)

  }

}

#' @return \code{NULL}
#'
#' @rdname seeker_snp_arch
#' @export
seeker_snp_arch.data.frame <- function(ID){
  message(paste(Sys.time(), 'Running `seeker_snp_arch` for data.frame'))
  ID1 <- as.matrix(ID)
  server <- "http://rest.ensembl.org/variation/human/"
  ligas <- paste0(server, ID1,"?pops=1;content-type=application/json")

  mydf <- ID[NULL,]
  for(i in 1:length(ligas)){

    r <- fromJSON(ligas[i])
    pop <- r[["mappings"]]
    pop_result <- data.frame(SNP = ID[i,], pop)

    mydf <- rbind(mydf, pop_result)
  }
  return(mydf)
}

