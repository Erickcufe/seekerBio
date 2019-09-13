#' seeker Single Nucleotide Polymorphism Context
#'
#' A generic function to search the type of variation & gene
#'
#' @param SNP A Single Nucleotide Polymorphism ID ("rs") in character or data.frame
#'
#' @return
#' A data.frame with Type of variation & gene of each SNP
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
  message(paste(Sys.time(), 'Running `seeker_snp_context` for factor'))

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
  message(paste(Sys.time(), 'Running `seeker_snp_context` for character'))

  if (length(SNP)==1){

    SNPs <- SNP
    SNPs_1 <- unlist(strsplit(SNPs, split = 'rs', fixed = T))
    SNPs_1 <- SNPs_1[SNPs_1 != ""]
    URL_dbSNP <- "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/"
    ligas <- paste0(URL_dbSNP,SNPs_1)
    vector_SNPs <- list()

    future::plan(multiprocess)

    contents <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                  .progress = TRUE)
    contents_1 <- purrr::transpose(contents)
    contents_request <- contents_1[["result"]]

    vector_SNPs <- list()
    # For context
    for (i in 1:length(contents_request)){

      request <- contents_request[[i]]$primary_snapshot_data$allele_annotations$assembly_annotation

      if (!is.null(request[[2]][["genes"]][[1]][["rnas"]][[1]][["protein"]][["sequence_ontology"]][[1]][["name"]])){

        context_SNP <- request[[2]][["genes"]][[1]][["rnas"]][[1]][["protein"]][["sequence_ontology"]][[1]][["name"]]

      } else {

        if(!is.null(request[[2]][["genes"]][[1]][["rnas"]][[1]][["sequence_ontology"]][[1]][["name"]])){

          context_SNP <- request[[2]][["genes"]][[1]][["rnas"]][[1]][["sequence_ontology"]][[1]][["name"]]

        } else {

          context_SNP <- request[[1]][["genes"]][[1]][["rnas"]][[1]][["sequence_ontology"]][[1]][["name"]]

        }
      }

      if(is.null(context_SNP)){

        context_SNP <- "UNKNOW"
      }

      vector_SNPs[[i]] <- context_SNP

    }

    SNPs_3 <- vector_SNPs[!sapply(vector_SNPs, is.null)]
    SNPs_3 <- data.frame(SNPs_3)
    SNPs_3 <- SNPs_3[1,]
    SNPs_3 <- t(SNPs_3)
    All_SNPs <- cbind(SNPs,as.character(SNPs_3))
    All_SNPs <- data.frame(All_SNPs)
    colnames(All_SNPs) <- c("SNP","CONTEXT")
    rownames(All_SNPs) <- NULL


    # For GENE

    vector_gene <- list()
    for (i in 1:length(contents_request)){

      request <- contents_request[[i]]$primary_snapshot_data$allele_annotations$assembly_annotation

      if (!is.null(request[[2]][["genes"]][[1]][["locus"]])){

        context_SNP <- request[[2]][["genes"]][[1]][["locus"]]

      } else {

        if(!is.null(request[[1]][["genes"]][[1]][["locus"]])){

          context_SNP <- request[[1]][["genes"]][[1]][["locus"]]

        } else {

          context_SNP <- "UNKNOW"

        }
      }


      vector_gene[[i]] <- context_SNP

    }

    SNPs_4 <- vector_gene[!sapply(vector_gene, is.null)]

    df_snp4 <- data.frame()
    for (i in 1:length(SNPs_4)) {

      if (length(SNPs_4[[i]])==1){

        tmp_df <- data.frame(CONTEXT=SNPs_4[[i]])
        df_snp4 <- rbind(df_snp4, tmp_df)

      } else {

        tmp_df <- data.frame(CONTEXT=paste(SNPs_4[[i]], collapse = ", "))
        df_snp4 <- rbind(df_snp4, tmp_df)
      }

    }

    genes <- cbind(SNPs_1,df_snp4)
    genes <- data.frame(genes)
    colnames(genes) <- c("SNP","GENE")
    rownames(genes) <- NULL
    context_gene <- cbind(All_SNPs, GENE=genes$GENE)
    rownames(context_gene) <- NULL

    return(context_gene)

  } else {

    df <- data.frame(SNP = SNP)
    pop_result <- seeker_snp_context(SNP=df)
    return(pop_result)


  }

}


#' @rdname seeker_snp_context
#' @export
seeker_snp_context.data.frame <- function(SNP){

  message(paste(Sys.time(), 'Running `seeker_snp_context` for data.frame'))


  SNPs <- as.matrix(SNP)

  SNPs_1 <- unlist(strsplit(SNPs, split = 'rs', fixed = T))
  SNPs_1 <- SNPs_1[SNPs_1 != ""]
  URL_dbSNP <- "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/"
  ligas <- paste0(URL_dbSNP,SNPs_1)

  future::plan(multiprocess)

  contents <- furrr::future_map(ligas, purrr::safely(jsonlite::fromJSON),
                                .progress = TRUE)

  contents_1 <- purrr::transpose(contents)
  contents_request <- contents_1[["result"]]

  vector_SNPs <- list()
  # For context
  for (i in 1:length(contents_request)){

    request <- contents_request[[i]]$primary_snapshot_data$allele_annotations$assembly_annotation

    if (!is.null(request[[2]][["genes"]][[1]][["rnas"]][[1]][["protein"]][["sequence_ontology"]][[1]][["name"]])){

      context_SNP <- request[[2]][["genes"]][[1]][["rnas"]][[1]][["protein"]][["sequence_ontology"]][[1]][["name"]]

    } else {

      if(!is.null(request[[2]][["genes"]][[1]][["rnas"]][[1]][["sequence_ontology"]][[1]][["name"]])){

        context_SNP <- request[[2]][["genes"]][[1]][["rnas"]][[1]][["sequence_ontology"]][[1]][["name"]]

      } else {

        context_SNP <- request[[1]][["genes"]][[1]][["rnas"]][[1]][["sequence_ontology"]][[1]][["name"]]

      }
    }

    if(is.null(context_SNP)){

      context_SNP <- "UNKNOW"
    }

    vector_SNPs[[i]] <- context_SNP

  }

  SNPs_3 <- vector_SNPs[!sapply(vector_SNPs, is.null)]
  SNPs_3 <- data.frame(SNPs_3)
  SNPs_3 <- SNPs_3[1,]
  SNPs_3 <- t(SNPs_3)
  All_SNPs <- cbind(SNPs,SNPs_3)
  All_SNPs <- data.frame(All_SNPs)
  colnames(All_SNPs) <- c("SNP","CONTEXT")
  rownames(All_SNPs) <- NULL


  # For GENE

  vector_gene <- list()
  for (i in 1:length(contents_request)){

    request <- contents_request[[i]]$primary_snapshot_data$allele_annotations$assembly_annotation

    if (!is.null(request[[2]][["genes"]][[1]][["locus"]])){

      context_SNP <- request[[2]][["genes"]][[1]][["locus"]]

    } else {

      if(!is.null(request[[1]][["genes"]][[1]][["locus"]])){

        context_SNP <- request[[1]][["genes"]][[1]][["locus"]]

      } else {

        context_SNP <- "UNKNOW"

      }
    }


    vector_gene[[i]] <- context_SNP

  }

  SNPs_4 <- vector_gene[!sapply(vector_gene, is.null)]

  df_snp4 <- data.frame()
  for (i in 1:length(SNPs_4)) {

    if (length(SNPs_4[[i]])==1){

      tmp_df <- data.frame(CONTEXT=SNPs_4[[i]])
      df_snp4 <- rbind(df_snp4, tmp_df)

    } else {

      tmp_df <- data.frame(CONTEXT=paste(SNPs_4[[i]], collapse = ", "))
      df_snp4 <- rbind(df_snp4, tmp_df)
    }

  }

  genes <- cbind(SNPs,df_snp4)
  genes <- data.frame(genes)
  colnames(genes) <- c("SNP","GENE")
  rownames(genes) <- NULL
  context_gene <- cbind(All_SNPs, GENE=genes$GENE)
  rownames(context_gene) <- NULL


  return(context_gene)
}

