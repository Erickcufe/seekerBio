seeker_snp_freq <- function(ID, study = "1000GENOMES:phase_3") {
  UseMethod("seeker_snp_freq")
}

seeker_snp_freq.character <- function(ID, study = "1000GENOMES:phase_3"){
  message(paste(Sys.time(), 'Running `seeker_snp_freq` for character'))

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

seeker_snp_freq.data.frame <- function(ID, study = "1000GENOMES:phase_3"){
  message(paste(Sys.time(), 'Running `seeker_snp_freq` for data.frame'))
  ID1 <- as.matrix(ID)
  server <- "http://rest.ensembl.org/variation/human/"
  ligas <- paste0(server, ID1,"?pops=1;content-type=application/json")

  mydf <- ID[NULL,]
  for(i in 1:length(ligas)){

    r <- fromJSON(ligas[i])
    pop <- r[["populations"]]

    seleccion <- stringr::str_detect(pop$population, study)
    SNP <- c(rep(ID, length(pop[seleccion,])))
    pop_result <- cbind(SNP = ID1[i], pop[seleccion,])
    pop_result$submission_id <- NULL

    mydf <- rbind(mydf, pop_result)
  }
  return(mydf)
}



