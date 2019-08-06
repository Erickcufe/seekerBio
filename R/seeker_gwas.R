seeker_gwas <- function(trait) {

  UseMethod("seeker_gwas")

}

seeker_gwas.character <- function(trait){
  # keyword_GWAS<-trait
  url_GWAS<-paste0("https://www.ebi.ac.uk/gwas/api/search/downloads?q=text:%22",
                       trait)
  GWAS<-read.delim(url_GWAS)
  separarcoma<-str_split(GWAS$REPORTED.GENE.S., pattern = fixed(","), n = Inf)
  separarcoma<-as.matrix(separarcoma)

  p_separar <- sapply(separarcoma, length)
  o_separar <- seq_len(max(p_separar))
  k_separar <- t(sapply(separarcoma, "[", i=o_separar))
  k_separar <- data.frame(k_separar, headers=F)
  names(k_separar)[1:length(k_separar)] <- c("GeneSymbol")


  binding<-data.frame(k_separar[,1])
  colnames(binding)<-"GeneSymbol"
  ENSEMBL_GWAS<-cbind(binding, GWAS)

  return(ENSEMBL_GWAS)
}


