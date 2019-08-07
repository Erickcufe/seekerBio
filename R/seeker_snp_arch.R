seeker_snp_arch <- function(ID){
  UseMethod("seeker_snp_arch")
}

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

