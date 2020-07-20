#' seeker_pathway_order
#'
#' This function order the pathways of reactome into categories
#'
#' @param pathways a data.frame product of the function seeker_gen_pathways()
#'
#' @return
#' A column with the classification pathway
#'
#' @source
#' https://reactome.org
#'
#' @author
#' Erick Cuevas Fernández
#'
#' Heriberto Manuel Rivera
#'
#' @examples
#'df <- data.frame(Gen=c("MMP12","MMP12","MMP12"), ID=c("R-HSA-1442490","R-HSA-1474228","R-HSA-1474244"), Path_name=c("Collagen degradation", "Degradation of the extracellular matrix", "Extracellular matrix organization"))
#'seeker_pathway_order(df)
#'
#' @rdname seeker_pathway_order
#' @export seeker_pathway_order
seeker_pathway_order <- function(pathways){
  # message(paste(Sys.time(), 'empezando a correr `seeker_pathway_order`'))


  #### LOOPS PARA AGRUPAR
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::metabolism_seeker)){
      pathways$General_pathway[i]<-"NA"
    }
  }

  #Loop para autophagy
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::autophagy_seeker)){
      if (pathways[i,2]==seekerBio::autophagy_seeker[j,1]){
        pathways$General_pathway[i]<-"Autophagy"
      }
    }
  }

  #Loop para metabolismo
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::metabolism_seeker)){
      if (pathways[i,2]==seekerBio::metabolism_seeker[j,1]){
        pathways$General_pathway[i]<-"Metabolism"
      }
    }
  }

  #Loop para signal transduction
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::signal_transduction_seeker)){
      if (pathways[i,2]==seekerBio::signal_transduction_seeker[j,1]){
        pathways$General_pathway[i]<-"Signal Transduction"
      }
    }
  }

  #Loop para Disease
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::disease_seeker)){
      if (pathways[i,2]==seekerBio::disease_seeker[j,1]){
        pathways$General_pathway[i]<-"Disease"
      }
    }
  }


  #Loop para Cell-Cell communication
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::cell_c_seeker)){
      if (pathways[i,2]==seekerBio::cell_c_seeker[j,1]){
        pathways$General_pathway[i]<-"Cell-Cell Communication"
      }
    }
  }

  #Loop para Cell Cycle
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::cell_cy_seeker)){
      if (pathways[i,2]==seekerBio::cell_cy_seeker[j,1]){
        pathways$General_pathway[i]<-"Cell Cycle"
      }
    }
  }

  #Loop para Programmed Cell death
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::cell_death_seeker)){
      if (pathways[i,2]==seekerBio::cell_death_seeker[j,1]){
        pathways$General_pathway[i]<-"Programmed Cell Death"
      }
    }
  }

  #Loop para Cellulaer responses to external stimuli
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::cell_re_seeker)){
      if (pathways[i,2]==seekerBio::cell_re_seeker[j,1]){
        pathways$General_pathway[i]<-"Cellular responses to externa stimuli"
      }
    }
  }

  #Loop para Chromatin organization
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::chromatin_seeker)){
      if (pathways[i,2]==seekerBio::chromatin_seeker[j,1]){
        pathways$General_pathway[i]<-"Chromatin organization"
      }
    }
  }

  #Loop para Circadian Clock
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::circadian_seeker)){
      if (pathways[i,2]==seekerBio::circadian_seeker[j,1]){
        pathways$General_pathway[i]<-"Circadian Clock"
      }
    }
  }

  #Loop para Developmental Biology
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::development_seeker)){
      if (pathways[i,2]==seekerBio::development_seeker[j,1]){
        pathways$General_pathway[i]<-"Developmental Biology"
      }
    }
  }

  #Loop para Digestion and absorption
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::digestion_seeker)){
      if (pathways[i,2]==seekerBio::digestion_seeker[j,1]){
        pathways$General_pathway[i]<-"Digestion and absorption"
      }
    }
  }

  #Loop para Extracellular matrix organization
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::extra_cell_matrix_seeker)){
      if (pathways[i,2]==seekerBio::extra_cell_matrix_seeker[j,1]){
        pathways$General_pathway[i]<-"Extracellular matrix organization"
      }
    }
  }

  #Loop para Gene expression (Transcription)
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::gene_expression_seeker)){
      if (pathways[i,2]==seekerBio::gene_expression_seeker[j,1]){
        pathways$General_pathway[i]<-"Gene expression (Transcription)"
      }
    }
  }

  #Loop para Hemostasis
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::hemostasis_seeker)){
      if (pathways[i,2]==seekerBio::hemostasis_seeker[j,1]){
        pathways$General_pathway[i]<-"Hemostasis"
      }
    }
  }

  #Loop para Immune System
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::immune_system_seeker)){
      if (pathways[i,2]==seekerBio::immune_system_seeker[j,1]){
        pathways$General_pathway[i]<-"Immune System"
      }
    }
  }

  #Loop para Metabolism of proteins
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::metabolism_proteins_seeker)){
      if (pathways[i,2]==seekerBio::metabolism_proteins_seeker[j,1]){
        pathways$General_pathway[i]<-"Metabolism of proteins"
      }
    }
  }

  #Loop para Metabolism of RNA
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::metabolism_rna_seeker)){
      if (pathways[i,2]==seekerBio::metabolism_rna_seeker[j,1]){
        pathways$General_pathway[i]<-"Metabolism of RNA"
      }
    }
  }

  #Loop para Mitophagy
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::mitophagy_seeker)){
      if (pathways[i,2]==seekerBio::mitophagy_seeker[j,1]){
        pathways$General_pathway[i]<-"Mitophagy"
      }
    }
  }

  #Loop para Muscle contraction
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::muscle_cons_seeker)){
      if (pathways[i,2]==seekerBio::muscle_cons_seeker[j,1]){
        pathways$General_pathway[i]<-"Muscle contraction"
      }
    }
  }

  #Loop para Neuronal System
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::neuronal_system_seeker)){
      if (pathways[i,2]==seekerBio::neuronal_system_seeker[j,1]){
        pathways$General_pathway[i]<-"Neuronal System"
      }
    }
  }

  #Loop para Organelle biogenesis and maintenance
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::organelle_biogenesis_seeker)){
      if (pathways[i,2]==seekerBio::organelle_biogenesis_seeker[j,1]){
        pathways$General_pathway[i]<-"Organelle biogenesis and maintenance"
      }
    }
  }

  #Loop para Protein Localization
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::protein_loc_seeker)){
      if (pathways[i,2]==seekerBio::protein_loc_seeker[j,1]){
        pathways$General_pathway[i]<-"Protein localization"
      }
    }
  }

  #Loop para DNA Repair
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::repairdna_seeker)){
      if (pathways[i,2]==seekerBio::repairdna_seeker[j,1]){
        pathways$General_pathway[i]<-"DNA Repair"
      }
    }
  }

  #Loop para DNA Replication
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::replication_seeker)){
      if (pathways[i,2]==seekerBio::replication_seeker[j,1]){
        pathways$General_pathway[i]<-"DNA Replication"
      }
    }
  }

  #Loop para Reproduction
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::reproduction_seeker)){
      if (pathways[i,2]==seekerBio::reproduction_seeker[j,1]){
        pathways$General_pathway[i]<-"Reproduction"
      }
    }
  }

  #Loop para Transport of small molecules
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::transport_small_molecules_seeker)){
      if (pathways[i,2]==seekerBio::transport_small_molecules_seeker[j,1]){
        pathways$General_pathway[i]<-"Transport of small molecules"
      }
    }
  }

  #Loop para Vesicle-mediated transport
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(seekerBio::vesicle_transport_seeker)){
      if (pathways[i,2]==seekerBio::vesicle_transport_seeker[j,1]){
        pathways$General_pathway[i]<-"Vesicle-mediated transport"
      }
    }
  }

    return(pathways)
}
