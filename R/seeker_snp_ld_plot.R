#' seeker_snp_ld_plot
#'
#' A function that plot the Linkage-Dsequilibrium information for a given Single Nucleotide Polymorphism
#'
#' @param SNP A Single Nucleotide Polymorphism ID ("rs") in character.
#' @param population_study Population for which to compute LD. Default "1000GENOMES:phase_3:MXL"
#' @param window_size Window size in kb. The maximum allowed value for the window size is 500 kb.
#' LD is computed for the given variant and all variants that are located within the specified window.
#' Default 500
#' @param d_prime Measure of LD. If D' is provided only return pairs of variants whose D' value is equal to or greater than the value provided.
#' Default 0
#' @param color_select Select a color for the plot; Default= "cyan"
#'
#' @importFrom
#' ggrepel geom_label_repel
#'
#' @import
#' ggplot2
#'
#' @author
#' Erick Cuevas Fernández
#'
#' Heriberto Manuel Rivera
#'
#' @return
#' A plot, in x axis the position in the chromosome, in y axis the r2 of linkage desequilibrium
#'
#'
#'
#' @examples
#' seeker_snp_ld_plot("rs7412", "1000GENOMES:phase_3:MXL", color_select = "green")
#'
#' seeker_snp_ld_plot("rs7412", "1000GENOMES:phase_3:SAS", color_select = "pink")
#'
#' @export seeker_snp_ld_plot
seeker_snp_ld_plot <- function(SNP, population_study="1000GENOMES:phase_3:MXL",
                               window_size = 500, d_prime = 0, color_select = "cyan") {

  LD <- seekerBio::seeker_snp_ld(SNP, population = population_study)

  if (length(LD)==0){
    stop("**Your SNP don´t have information in https://rest.ensembl.org, sorry**", call. = FALSE)
  }
  variations_arch <- seekerBio::seeker_snp_arch(LD$variation2)
  variation_snp <- seekerBio::seeker_snp_arch(SNP)

  to_plot <- merge.data.frame(LD, variations_arch, by.x = "variation2",
                              by.y = "SNP")

  variation_snp$r2 <- 1

  label1 <- paste("Chr", variation_snp$seq_region_name)

  ggplot(data = to_plot, aes(x = start, y = as.numeric(r2), colour = as.numeric(d_prime))) +
    geom_point(alpha = 0.75, size = 2) + geom_point(data = variation_snp,
                                                    aes(x = start, y = as.numeric(r2)),
                                                    size= 4,
                                                    color= "red",
                                                    alpha= 0.8) +
    theme_minimal() + geom_area(aes(start), data=to_plot, alpha= 0.2,
                                fill= color_select, color="black") +
    scale_color_gradient(low="green", high="red") +
    ggrepel::geom_label_repel(aes(label=ifelse(r2>0.6,as.character(variation2),'')),box.padding= 0.2,
                              point.padding = 0.3,
                              segment.color = 'grey50') +
    theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size = 20),
          axis.text.y =element_text(size=20),axis.text.y.right = element_text(size=20),
          plot.title = element_text(hjust = 0.5, size=30),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.title.x=element_text(size=20),
          axis.title.y = element_text(size=20)) +
    xlab (label1) + ylab("r2") +
    geom_label(data= variation_snp,
               aes(x = start, y= as.numeric(r2)), label= variation_snp$SNP,
               nudge_x = 0.8, nudge_y = 0.02,
               color="red") +
    ggtitle(population_study) + labs(colour = "LD")

}
