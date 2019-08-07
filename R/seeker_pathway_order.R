#' seeker_pathway_order
#'
#' This function order the pathways of reactome into categories
#'
#' @param pathways a data.frame product of the function seeker_gen_pathways()
#'
#' @return
#' A column with the classification pathway
#'
#'
#' @examples
#'df <- data.frame(Gen=c("MMP12","MMP12","MMP12"), ID=c("R-HSA-1442490","R-HSA-1474228","R-HSA-1474244"), Path_name=c("Collagen degradation", "Degradation of the extracellular matrix", "Extracellular matrix organization"))
#'seeker_pathway_order(df)
#'
#' @rdname seeker_pathway_order
#' @export seeker_pathway_order
seeker_pathway_order <- function(pathways){
  message(paste(Sys.time(), 'empezando a correr `seeker_pathway_order`'))
  #### LOOPS PARA AGRUPAR
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(metabolism_seeker)){
      pathways$General_pathway[i]<-"NA"
    }
  }

  #Loop para autophagy
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(autophagy_seeker)){
      if (pathways[i,2]==autophagy_seeker[j,1]){
        pathways$General_pathway[i]<-"Autophagy"
      }
    }
  }

  #Loop para metabolismo
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(metabolism_seeker)){
      if (pathways[i,2]==metabolism_seeker[j,1]){
        pathways$General_pathway[i]<-"Metabolism"
      }
    }
  }

  #Loop para signal transduction
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(signal_transduction_seeker)){
      if (pathways[i,2]==signal_transduction_seeker[j,1]){
        pathways$General_pathway[i]<-"Signal Transduction"
      }
    }
  }

  #Loop para Disease
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(disease_seeker)){
      if (pathways[i,2]==disease_seeker[j,1]){
        pathways$General_pathway[i]<-"Disease"
      }
    }
  }


  #Loop para Cell-Cell communication
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(cell_c_seeker)){
      if (pathways[i,2]==cell_c_seeker[j,1]){
        pathways$General_pathway[i]<-"Cell-Cell Communication"
      }
    }
  }

  #Loop para Cell Cycle
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(cell_cy_seeker)){
      if (pathways[i,2]==cell_cy_seeker[j,1]){
        pathways$General_pathway[i]<-"Cell Cycle"
      }
    }
  }

  #Loop para Programmed Cell death
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(cell_death_seeker)){
      if (pathways[i,2]==cell_death_seeker[j,1]){
        pathways$General_pathway[i]<-"Programmed Cell Death"
      }
    }
  }

  #Loop para Cellulaer responses to external stimuli
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(cell_re_seeker)){
      if (pathways[i,2]==cell_re_seeker[j,1]){
        pathways$General_pathway[i]<-"Cellular responses to externa stimuli"
      }
    }
  }

  #Loop para Chromatin organization
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(chromatin_seeker)){
      if (pathways[i,2]==chromatin_seeker[j,1]){
        pathways$General_pathway[i]<-"Chromatin organization"
      }
    }
  }

  #Loop para Circadian Clock
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(circadian_seeker)){
      if (pathways[i,2]==circadian_seeker[j,1]){
        pathways$General_pathway[i]<-"Circadian Clock"
      }
    }
  }

  #Loop para Developmental Biology
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(development_seeker)){
      if (pathways[i,2]==development_seeker[j,1]){
        pathways$General_pathway[i]<-"Developmental Biology"
      }
    }
  }

  #Loop para Digestion and absorption
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(digestion_seeker)){
      if (pathways[i,2]==digestion_seeker[j,1]){
        pathways$General_pathway[i]<-"Digestion and absorption"
      }
    }
  }

  #Loop para Extracellular matrix organization
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(extra_cell_matrix_seeker)){
      if (pathways[i,2]==extra_cell_matrix_seeker[j,1]){
        pathways$General_pathway[i]<-"Extracellular matrix organization"
      }
    }
  }

  #Loop para Gene expression (Transcription)
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(gene_expression_seeker)){
      if (pathways[i,2]==gene_expression_seeker[j,1]){
        pathways$General_pathway[i]<-"Gene expression (Transcription)"
      }
    }
  }

  #Loop para Hemostasis
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(hemostasis_seeker)){
      if (pathways[i,2]==hemostasis_seeker[j,1]){
        pathways$General_pathway[i]<-"Hemostasis"
      }
    }
  }

  #Loop para Immune System
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(immune_system_seeker)){
      if (pathways[i,2]==immune_system_seeker[j,1]){
        pathways$General_pathway[i]<-"Immune System"
      }
    }
  }

  #Loop para Metabolism of proteins
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(metabolism_proteins_seeker)){
      if (pathways[i,2]==metabolism_proteins_seeker[j,1]){
        pathways$General_pathway[i]<-"Metabolism of proteins"
      }
    }
  }

  #Loop para Metabolism of RNA
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(metabolism_rna_seeker)){
      if (pathways[i,2]==metabolism_rna_seeker[j,1]){
        pathways$General_pathway[i]<-"Metabolism of RNA"
      }
    }
  }

  #Loop para Mitophagy
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(mitophagy_seeker)){
      if (pathways[i,2]==mitophagy_seeker[j,1]){
        pathways$General_pathway[i]<-"Mitophagy"
      }
    }
  }

  #Loop para Muscle contraction
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(muscle_cons_seeker)){
      if (pathways[i,2]==muscle_cons_seeker[j,1]){
        pathways$General_pathway[i]<-"Muscle contraction"
      }
    }
  }

  #Loop para Neuronal System
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(neuronal_system_seeker)){
      if (pathways[i,2]==neuronal_system_seeker[j,1]){
        pathways$General_pathway[i]<-"Neuronal System"
      }
    }
  }

  #Loop para Organelle biogenesis and maintenance
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(organelle_biogenesis_seeker)){
      if (pathways[i,2]==organelle_biogenesis_seeker[j,1]){
        pathways$General_pathway[i]<-"Organelle biogenesis and maintenance"
      }
    }
  }

  #Loop para Protein Localization
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(protein_loc_seeker)){
      if (pathways[i,2]==protein_loc_seeker[j,1]){
        pathways$General_pathway[i]<-"Protein localization"
      }
    }
  }

  #Loop para DNA Repair
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(repairdna_seeker)){
      if (pathways[i,2]==repairdna_seeker[j,1]){
        pathways$General_pathway[i]<-"DNA Repair"
      }
    }
  }

  #Loop para DNA Replication
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(replication_seeker)){
      if (pathways[i,2]==replication_seeker[j,1]){
        pathways$General_pathway[i]<-"DNA Replication"
      }
    }
  }

  #Loop para Reproduction
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(reproduction_seeker)){
      if (pathways[i,2]==reproduction_seeker[j,1]){
        pathways$General_pathway[i]<-"Reproduction"
      }
    }
  }

  #Loop para Transport of small molecules
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(transport_small_molecules_seeker)){
      if (pathways[i,2]==transport_small_molecules_seeker[j,1]){
        pathways$General_pathway[i]<-"Transport of small molecules"
      }
    }
  }

  #Loop para Vesicle-mediated transport
  for(i in 1:nrow(pathways)){
    for (j in 1:nrow(vesicle_transport_seeker)){
      if (pathways[i,2]==vesicle_transport_seeker[j,1]){
        pathways$General_pathway[i]<-"Vesicle-mediated transport"
      }
    }
  }

    return(pathways)
}
