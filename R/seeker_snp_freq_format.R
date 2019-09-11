#' seeker_snp_freq_format
#'
#' @param data A data.frame, output of seeker_snp_freq()
#'
#' @return
#' A data.frame with the allel frequency in a horizontal format, with sorted frequencies.
#'
#' @author
#' Erick Cuevas-Fernández
#'
#'
#' @examples
#'
#' df <- seeker_snp_freq("rs56116432")
#'
#' seeker_snp_freq_format(df)
#'
#' @rdname seeker_snp_freq_format
#' @export seeker_snp_freq_format
seeker_snp_freq_format <- function(data){

  SNPs <- unique(data$SNP)
  list_snp <- list()

  for (i in 1:length(SNPs)){

    df <- data.frame(data[data$SNP %in% SNPs[i],])
    df$allele_count <- NULL

    list_snp[[i]] <- df

  }


  # Start to classified
  completos <- list()
  tres_alelos <- list()
  incompletos <- list()

  for (i in 1:length(list_snp)){

    df <- list_snp[[i]]

    if (nrow(df) == 64){

      completos[[i]] <- list_snp[[i]]

    }

    if (nrow(df) < 64){

      incompletos[[i]] <- list_snp[[i]]

    }

    if(nrow(df) > 64) {

      tres_alelos[[i]] <- list_snp[[i]]

    }

  }



  completos[sapply(completos, is.null)] <- NULL

  if(length(tres_alelos)!=0){

    tres_alelos[sapply(tres_alelos, is.null)] <- NULL

  }

  incompletos[sapply(incompletos, is.null)] <- NULL

  print("estamos en este paso")
  #Arrange Completos, incompletos and tres_alelos
  print(length(completos))

  names_pop <- c("Global", "Global_MAF", "ACB", "ACB_MAF", "AFR", "AFR_MAF",
                 "AMR", "AMR_MAF", "ASW", "ASW_MAF", "BEB", "BEB_MAF", "CDX",
                 "CDX_MAF", "CEU", "CEU_MAF", "CHB", "CHB_MAF", "CHS", "CHS_MAF",
                 "CLM", "CLM_MAF", "EAS", "EAS_MAF", "ESN", "ESN_MAF", "EUR",
                 "EUR_MAF", "FIN", "FIN_MAF", "GBR", "GBR_MAF", "GIH", "GIH_MAF",
                 "IBS", "IBS_MAF", "ITU", "ITU_MAF", "JPT", "JPT_MAF", "KHV",
                 "KHV_MAF", "LWK", "LWK_MAF", "GWD", "GWD_MAF", "MSL", "MSL_MAF",
                 "MXL", "MXL_MAF", "PEL", "PEL_MAF", "PJL", "PJL_MAF", "PUR",
                 "PUR_MAF", "SAS", "SAS_MAF", "STU", "STU_MAF", "TSI", "TSI_MAF",
                 "YRI", "YRI_MAF")

  # COMPLETE CASES
  mydf_complete <- data.frame()

  for (i in 1:length(completos)){


    df_completos <- completos[[i]]

    for(j in seq(from= 1, to=nrow(df_completos), by=2)){

      ordered_freq <- sort(df_completos[c(j,j+1), 4], decreasing = TRUE)
      df_completos[c(j,j+1), 4] <- ordered_freq

    }

    df_completos_1 <- data.frame(t(df_completos[,c(3,4)]))
    colnames(df_completos_1) <- names_pop
    df_completos_1 <- df_completos_1[-1,]

    df_completos_2 <- cbind(SNP = df_completos[1,1], df_completos_1,
                            Ancestral = df_completos[1,2],
                            Minor = df_completos[2,2])

    mydf_complete <- rbind(mydf_complete, df_completos_2)


  }


  # INCOMPLETE CASES

  for (i in 1:length(incompletos)){

    incomplete_df <- incompletos[[i]]

    mydf_incomplete <- data.frame()
    for (j in seq(from=1, to= nrow(incomplete_df), by=2)){

      if(j == nrow(incomplete_df)){
        break
      }


      if (incomplete_df[c(j),3]==incomplete_df[c(j+1),3]){

        mydf_incomplete <- rbind(mydf_incomplete,
                                 incomplete_df[c(j, j+1),] )

      } else{

        if(incomplete_df[j,4]< 1){

          tmp_inc <- data.frame(incomplete_df[j,c(1,2,3)],
                                frequency=(1-incomplete_df[j,4]))

          mydf_incomplete <- rbind(mydf_incomplete,
                                   incomplete_df[j,],
                                   tmp_inc)

          next

        }

        tmp_inc <- data.frame(incomplete_df[j,c(1,2,3)],
                              frequency=0)
        mydf_incomplete <- rbind(mydf_incomplete,
                                 incomplete_df[j,],
                                 tmp_inc)

      }
    }

    incompletos[[i]] <- mydf_incomplete

  }

  print(incompletos)
  print("casi terminamos")


  example_pop <- completos[[1]]

  for(i in 1:length(incompletos)){

    dummy_data <- incompletos[[i]]
    missing_pop <- example_pop$population %in% dummy_data$population

    tmp_pop <- data.frame()
    for(j in seq(from=1, to=64, by=2)){

      if(missing_pop[j]==TRUE){

        pop_wa <- example_pop[j,3]
        pop_catch <- dummy_data[dummy_data$population %in% pop_wa, c(1,2,3,4)]
        tmp_pop <- rbind(tmp_pop, pop_catch)

      } else {

        pop_catch <- data.frame(SNP= dummy_data[c(1,2),1],
                                allele= dummy_data[c(1,2),2],
                                population = example_pop[c(j,j+1), 3],
                                frequency=c(1,0))
        tmp_pop <- rbind(tmp_pop, pop_catch)

      }

    }

    incompletos[[i]] <- tmp_pop
  }

  # return(incompletos)


  mydf_non <- data.frame()

  for (i in 1:length(incompletos)){

    df_incompletos <- incompletos[[i]]

    for(j in seq(from= 1, to=nrow(df_incompletos), by=2)){

      ordered_freq <- sort(df_incompletos[c(j,j+1), 4], decreasing = TRUE)
      df_incompletos[c(j,j+1), 4] <- ordered_freq

    }

    df_completos_1 <- data.frame(t(df_incompletos[,c(3,4)]))
    colnames(df_completos_1) <- names_pop
    df_completos_1 <- df_completos_1[-1,]

    df_completos_2 <- cbind(SNP = df_completos[1,1], df_completos_1,
                            Ancestral = df_completos[1,2],
                            Minor = df_completos[2,2])

    mydf_non <- rbind(mydf_non, df_completos_2)


  }

  mydf_all <- rbind(mydf_complete, mydf_non)

  return(mydf_all)



}
