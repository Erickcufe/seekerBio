#' seeker_haploview
#'
#' A function that plot a correlation r2 between all SNPs with a Linkage Desequilibrium upper 0.5
#'
#' @param SNP A Single Nucleotide Polymorphism ID ("rs") in character.
#' @param population_study Population for which to compute LD. Default "1000GENOMES:phase_3:MXL"
#'
#'
#'
#' @import
#' ggplot2
#'
#' @import
#' dplyr
#'
#' @author
#' Erick Cuevas Fernández
#'
#' Heriberto Manuel Rivera
#'
#' @return
#' A correlation r2 plot between all SNPs with a Linkage Desequilibrium upper 0.5
#'
#'
#'
#' @examples
#' seeker_haploview("rs7412")
#'
#' seeker_haploview("rs7412", population_study = "1000GENOMES:phase_3:FIN")
#'
#' @export seeker_haploview
seeker_haploview <- function(SNP, population_study="1000GENOMES:phase_3:MXL",
                             window_size = 500, d_prime = 0, color_select = "cyan"){


  LD <- seekerBio::seeker_snp_ld(SNP, population = population_study)

  if (length(LD)==0){
    stop("**Your SNP don´t have information in https://rest.ensembl.org, sorry**", call. = FALSE)
  }

  variation_snp <- seekerBio::seeker_snp_arch(SNP)

  most_LD <- LD %>% dplyr::filter(r2 >= 0.6)

  LD_correlation <- seekerBio::seeker_snp_ld(most_LD$variation2,
                                             population = population_study)

  finder <- LD_correlation[LD_correlation$variation2 %in% most_LD$variation2, ]

  finder_1 <- rbind(most_LD, finder)

  finder_2 <- data.frame(v1= finder_1$variation1, v2= finder_1$variation2,
                         r2= finder_1$r2)
  tmp_sn1 <- data.frame(v1=c(SNP, most_LD$variation2),
                        v2=c(SNP, most_LD$variation2),
                        r2 = c(1))
  tmp_sn2 <- data.frame(v1 = most_LD$variation2,
                        v2 = most_LD$variation1,
                        r2 = most_LD$r2)


  finder_3 <- rbind(tmp_sn1, finder_2, tmp_sn2)

  finder_3$r2 <- as.character(finder_3$r2)
  finder_3$r2 <- as.numeric(finder_3$r2)

  label1 <- paste("Chr", variation_snp$seq_region_name)

  a <- ifelse(levels(finder_3$v2) == SNP, "red", "black")

  ggplot2::ggplot(finder_3, aes(x = v1, y= v2)) +
    geom_tile(aes(fill=r2), color= 'white') +
    scale_fill_gradient2(low="blue", high="red", mid="yellow",
                                        midpoint = 0.5,
                                        limit=c(0.1,1),
                                        space = "Lab",
                                        name="r2") +
    theme_minimal()  +
    coord_fixed() + theme(panel.grid = element_blank()) +
    theme(axis.text.x = element_text(angle = 60, hjust = 0.9, size = 20, color = a),
          axis.text.y =element_text(size=20, color=a),
          axis.text.y.right = element_text(size=20),
          axis.ticks=element_blank(),
          axis.line=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_line(color='#eeeeee'),
          plot.title = element_text(hjust = 0.5, size=30),
          legend.title = element_text(size = 20),
          axis.title.x=element_text(size=20),
          axis.title.y=element_text(size=20)) +
    ggtitle(paste(label1, ":", population_study)) + xlab ("SNPs") + ylab("SNPs")

}
