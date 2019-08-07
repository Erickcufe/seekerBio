#' seeker gwas
#'
#'seeker_gwas is a generic function that allow to download a data.frame of all Genome-wide
#'Association Study of the input trait
#'
#' @param trait A keyword to search in GWAS
#'
#' @return
#' A data.frame with all the information from https://www.ebi.ac.uk/gwas/home
#'
#' @importFrom
#' stringr str_split
#'
#' @examples
#' seeker_gwas("multiple sclerosis")
#' seeker_gwas("fontotemporal")
#' @rdname seeker_gwas
#' @export seeker_gwas
seeker_gwas <- function(trait) {

  UseMethod("seeker_gwas")

}

#' @return \code{NULL}
#'
#' @rdname seeker_gwas
#' @export
seeker_gwas.character <- function(trait){
  message(paste(Sys.time(), 'empezando a correr `seeker_gen_pathway`'))
  trait <-  strsplit(trait, " ")

  if (!is.na(trait[[1]][3])){
    stop("You must specified your trait with TWO words")

    } else {

    if (!is.na(trait[[1]][2])) {

      print("two words")
      # trait <- "multiple sclerosis"

      keyword_GWAS <- trait[[1]][1]
      keyword_2GWAS <- trait[[1]][2]

      url_GWAS <- paste0("https://www.ebi.ac.uk/gwas/api/search/downloads?q=text:%22",
                       keyword_GWAS,"%20",keyword_2GWAS,
                       "%22&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&traitfilter[]=&genotypingfilter[]=&dateaddedfilter=&efo=true&facet=association")

      GWAS <- read.delim(url_GWAS)
      separarcoma <- stringr::str_split(GWAS$REPORTED.GENE.S., pattern = stringr::fixed(","), n = Inf)
      separarcoma <- as.matrix(separarcoma)

      p_separar <- sapply(separarcoma, length)
      o_separar <- seq_len(max(p_separar))
      k_separar <- t(sapply(separarcoma, "[", i=o_separar))
      k_separar <- data.frame(k_separar, headers=F)
      names(k_separar)[1:length(k_separar)] <- c("GeneSymbol")


      binding <- data.frame(k_separar[,1])

      separarcoma_2 <- stringr::str_split(binding$k_separar...1., pattern = stringr::fixed(" "), n = Inf)
      separarcoma_2 <- as.matrix(separarcoma_2)

      p_separar_2 <- sapply(separarcoma_2, length)
      o_separar_2 <- seq_len(max(p_separar_2))
      k_separar_2 <- t(sapply(separarcoma_2, "[", i=o_separar))
      k_separar_2 <- data.frame(k_separar_2, headers=F)
      names(k_separar_2)[1:length(k_separar_2)] <- c("GeneSymbol")
      binding_2 <- data.frame(k_separar_2[,1])

      colnames(binding_2) <- "GeneSymbol"
      ENSEMBL_GWAS <- cbind(binding_2, GWAS)


    } else {

      if (is.na(trait[[1]][2])){

        print("one_word")
        url_GWAS <- paste0("https://www.ebi.ac.uk/gwas/api/search/downloads?q=text:%22",
                         trait,
                         "%22&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&traitfilter[]=&genotypingfilter[]=&dateaddedfilter=&efo=true&facet=association")
        GWAS <- read.delim(url_GWAS)
        separarcoma <- stringr::str_split(GWAS$REPORTED.GENE.S., pattern = stringr::fixed(","), n = Inf)
        separarcoma <- as.matrix(separarcoma)

        p_separar <- sapply(separarcoma, length)
        o_separar <- seq_len(max(p_separar))
        k_separar <- t(sapply(separarcoma, "[", i=o_separar))
        k_separar <- data.frame(k_separar, headers=F)
        names(k_separar)[1:length(k_separar)] <- c("GeneSymbol")


        binding <- data.frame(k_separar[,1])


        separarcoma_2 <- stringr::str_split(binding$k_separar...1., pattern = stringr::fixed(" "), n = Inf)
        separarcoma_2 <- as.matrix(separarcoma_2)

        p_separar_2 <- sapply(separarcoma_2, length)
        o_separar_2 <- seq_len(max(p_separar_2))
        k_separar_2 <- t(sapply(separarcoma_2, "[", i=o_separar))
        k_separar_2 <- data.frame(k_separar_2, headers=F)
        names(k_separar_2)[1:length(k_separar_2)] <- c("GeneSymbol")
        binding_2 <- data.frame(k_separar_2[,1])

        colnames(binding_2) <- "GeneSymbol"
        ENSEMBL_GWAS <- cbind(binding_2, GWAS)

      }
    }
  }

  return(ENSEMBL_GWAS)
}


