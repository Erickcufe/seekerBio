#' seeker_snp_freq_format
#'
#' @param data A data.frame, output of seeker_snp_freq()
#'
#' @return
#' A list with the allele frequency in a horizontal format, with sorted frequencies. Two loci and three loci
#'
#' @author
#' Erick Cuevas-Fernández
#'
#' @import
#' dplyr
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


  if(length(completos)!=0){

  completos[sapply(completos, is.null)] <- NULL

  }

  if(length(tres_alelos)!=0){

    tres_alelos[sapply(tres_alelos, is.null)] <- NULL

  }


  if(length(incompletos)!=0){

    incompletos[sapply(incompletos, is.null)] <- NULL

  }
  # print("estamos en este paso")
  #Arrange Completos, incompletos and tres_alelos
  # print(completos[[1]])

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
  mydf_complete <- data.frame()
  # COMPLETE CASES
  if(length(completos)!=0){
  mydf_complete <- data.frame()

  for (i in 1:length(completos)){


    df_completos <- completos[[i]]

    for(j in seq(from= 1, to=nrow(df_completos), by=2)){

      ordered_freq <- sort(df_completos$frequency[c(j,j+1)], decreasing = TRUE)
      df_completos$frequency[c(j,j+1)] <- ordered_freq

    }

    tmp_completos <- data.frame(population = df_completos$population,
                                frequency = df_completos$frequency)
    df_completos_1 <- data.frame(t(tmp_completos))
    colnames(df_completos_1) <- names_pop
    df_completos_1 <- df_completos_1[-1,]

    df_completos_2 <- cbind(SNP = df_completos$SNP[1], df_completos_1,
                            Ancestral = df_completos$allele[1],
                            Minor = df_completos$allele[2])

    mydf_complete <- rbind(mydf_complete, df_completos_2)

    }
  }


  # INCOMPLETE CASES
  if(length(incompletos)!=0){
    for (i in 1:length(incompletos)){

      incomplete_df <- incompletos[[i]]
      menor <- incomplete_df[which(incomplete_df$frequency < 0.3)[1], "allele"]
      mayor <- incomplete_df[which(incomplete_df$frequency == 1)[1], "allele"]
      mydf_incomplete <- data.frame()
      for (j in seq(from=1, to= nrow(incomplete_df), by=2)){

        if(j == nrow(incomplete_df)){
          break
        }


        if (incomplete_df$population[j]==incomplete_df$population[c(j+1)]){

          tmp_incomplete <- data.frame(SNP=incomplete_df$SNP[c(j, j+1)],
                                       allele=incomplete_df$allele[c(j, j+1)],
                                       population=incomplete_df$population[c(j, j+1)],
                                       frequency=incomplete_df$frequency[c(j, j+1)])
          mydf_incomplete <- rbind(mydf_incomplete,
                                   tmp_incomplete)

        } else {

          if(incomplete_df$frequency[j]< 1){

            tmp_incomplete <- data.frame(SNP=incomplete_df$SNP[j],
                                         allele=mayor,
                                         population=incomplete_df$population[j],
                                         frequency=incomplete_df$frequency[j])

            tmp_inc <- data.frame(SNP=incomplete_df$SNP[j],
                                  allele=menor,
                                  population=incomplete_df$population[j],
                                  frequency=(1-incomplete_df$frequency[j]))

            mydf_incomplete <- rbind(mydf_incomplete,
                                     tmp_incomplete,
                                     tmp_inc)

            next

          }

          tmp_incomplete <- data.frame(SNP=incomplete_df$SNP[j],
                                       allele=mayor,
                                       population=incomplete_df$population[j],
                                       frequency=incomplete_df$frequency[j])

          tmp_inc <- data.frame(SNP=incomplete_df$SNP[j],
                                allele=menor,
                                population=incomplete_df$population[j],
                                frequency=0)
          mydf_incomplete <- rbind(mydf_incomplete,
                                   tmp_incomplete,
                                   tmp_inc)

        }
      }

      incompletos[[i]] <- mydf_incomplete

    }
  }


  example_pop <- seekerBio::populations_1000G
  if(length(incompletos)!=0){
    for(i in 1:length(incompletos)){

      dummy_data <- incompletos[[i]]
      missing_pop <- example_pop %in% dummy_data$population

      tmp_pop <- data.frame()
      for(j in seq(from=1, to=64, by=2)){

        if(missing_pop[j]==TRUE){

          pop_wa <- example_pop[j]
          pop_catch <- dummy_data[dummy_data$population %in% pop_wa, c(1,2,3,4)]
          tmp_pop <- rbind(tmp_pop, pop_catch)

        } else {

          pop_catch <- data.frame(SNP= dummy_data[c(1,2),1],
                                  allele= dummy_data[c(1,2),2],
                                  population = example_pop[c(j,j+1)],
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

      if(nrow(incompletos[[i]]) > 64){
        next
      }

      for(j in seq(from= 1, to=nrow(df_incompletos), by=2)){

        ordered_freq <- sort(df_incompletos[c(j,j+1), 4], decreasing = TRUE)
        df_incompletos[c(j,j+1), 4] <- ordered_freq

      }

      df_completos_1 <- data.frame(t(df_incompletos[,c(3,4)]))
      colnames(df_completos_1) <- names_pop
      df_completos_1 <- df_completos_1[-1,]

      df_completos_2 <- cbind(SNP = df_incompletos[1,1], df_completos_1,
                              Ancestral = df_incompletos$allele[5],
                              Minor = df_incompletos$allele[6])

      mydf_non <- rbind(mydf_non, df_completos_2)


    }
  }

  if(length(incompletos)!=0){
    mydf_all <- rbind(mydf_complete, mydf_non)
  } else {
    mydf_all <- mydf_complete

  }

  if(length(incompletos)!=0){
  for (i in 2:65) {

    mydf_all[,i] <- as.numeric(as.character(mydf_all[,i]))
    }
  }

  rownames(mydf_all) <- NULL

  if(length(tres_alelos)>1){

    data_for_three <- data.frame()
    bb <- data.frame()
    bb_3 <- data.frame()
    flag <- "1000GENOMES:phase_3:ALL"

    for (i in 1:length(tres_alelos)) {

      if(sum(tres_alelos[[i]]$population == flag) == 3){
        pp_3 <- tres_alelos[[i]]
        bb_3 <- rbind(bb_3, pp_3)
        next
      }

      pp <- tres_alelos[[i]]
      pp$population <- as.factor(pp$population)
      stoped <- "1000GENOMES:phase_3:YRI"
      indice <- which(pp$population==stoped)
      if(length(indice) == 2){
        if((indice[2] - indice[1]) > 1){
          pp <- pp[1:indice[1], ]
        } else {
          if((indice[2] - indice[1]) == 1){
            pp <- pp[1:indice[2], ]
          }
        }
      }

      if(length(indice) > 2){
        if((indice[2] - indice[1]) == 1){
          pp <- pp[1:indice[2], ]
        } else{
          pp <- pp[1:indice[1], ]
        }
      }
      bb <- rbind(bb, pp)
    }
    if(length(bb) > 0){
      data_for_three <- seekerBio::seeker_snp_freq_format(bb)
      mydf_all <- rbind(mydf_all, data_for_three[[1]])
    }

    mydf_all <- list(Two_loci = mydf_all, Three_loci = bb_3)
  } else {
    mydf_all <- list(Two_loci = mydf_all)
  }
  # message(paste("You have",length(tres_alelos), "SNPs with 3 allel frequency"))
  return(mydf_all)


}
