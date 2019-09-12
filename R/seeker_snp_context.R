seeker_snp_context <- function(){
  SNPs <- as.character(alz_ordered$SNP)

  SNPs_1 <- unlist(strsplit(SNPs, split = 'rs', fixed = T))
  SNPs_1 <- SNPs_1[SNPs_1 != ""]
  #SNPs_1<-SNPs_1[-1]

  URL_dbSNP <- "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/"
  vector_SNPs <- list()
  library(jsonlite)

  for (i in 1:length(SNPs_1)) {
    URL_paste <- paste0(URL_dbSNP,SNPs_1[1])
    URL_search <- fromJSON(URL_paste)
    URL_search <- URL_search$primary_snapshot_data$allele_annotations$assembly_annotation
    if (!is.null(URL_search[[2]][["genes"]][[1]][["rnas"]][[1]][["protein"]][["sequence_ontology"]][[1]][["name"]])){
      data_SNP <- URL_search[[2]][["genes"]][[1]][["rnas"]][[1]][["protein"]][["sequence_ontology"]][[1]][["name"]]

    } else {
      data_SNP <- URL_search[[1]][["genes"]][[1]][["rnas"]][[1]][["sequence_ontology"]][[1]][["name"]]

    }
    vector_SNPs[[i]] <- data_SNP
    print(i)

  }


  selector<-vector()
  for (i in 1:length(vector_SNPs)){
    selector[i] <- !is.null(vector_SNPs[[i]])
  }

  SNPs <- SNPs[-1,]
  SNPS_2 <- SNPs[selector]
  SNPs_3 <- vector_SNPs[!sapply(vector_SNPs, is.null)]
  SNPs_3 <- data.frame(SNPs_3)
  SNPs_3 <- SNPs_3[1,]
  SNPs_3 <- t(SNPs_3)
  All_SNPs <- cbind(SNPS_2,SNPs_3)
  All_SNPs <- data.frame(All_SNPs)

  colnames(All_SNPs) <- c("SNPS","CONTEXT")

}

