

#### LOOPS PARA AGRUPAR
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(metabolism)){
    Reactome_pathways1$Clasificador[i]<-"NA"
  }
}

#Loop para metabolismo
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(metabolism)){
    if (Reactome_pathways1[i,2]==metabolism[i,1]){
      Reactome_pathways1$Clasificador[i]<-"Metabolism"
    }
  }
}

#Loop para signal transduction
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(signal_transduction)){
    if (Reactome_pathways1[i,2]==signal_transduction[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Signal Transduction"
    }
  }
}

#Loop para Disease
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(disease)){
    if (Reactome_pathways1[i,2]==disease[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Disease"
    }
  }
}


#Loop para Cell-Cell communication
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(cell_c)){
    if (Reactome_pathways1[i,2]==cell_c[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Cell-Cell Communication"
    }
  }
}

#Loop para Cell Cycle
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(cell_cy)){
    if (Reactome_pathways1[i,2]==cell_cy[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Cell Cycle"
    }
  }
}

#Loop para Programmed Cell death
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(cell_death)){
    if (Reactome_pathways1[i,2]==cell_death[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Programmed Cell Death"
    }
  }
}

#Loop para Cellulaer responses to external stimuli
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(cell_re)){
    if (Reactome_pathways1[i,2]==cell_re[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Cellular responses to externa stimuli"
    }
  }
}

#Loop para Chromatin organization
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(chromatin)){
    if (Reactome_pathways1[i,2]==chromatin[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Chromatin organization"
    }
  }
}

#Loop para Circadian Clock
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(circadian)){
    if (Reactome_pathways1[i,2]==circadian[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Circadian Clock"
    }
  }
}

#Loop para Developmental Biology
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(development)){
    if (Reactome_pathways1[i,2]==development[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Developmental Biology"
    }
  }
}

#Loop para Digestion and absorption
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(digestion)){
    if (Reactome_pathways1[i,2]==digestion[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Digestion and absorption"
    }
  }
}

#Loop para Extracellular matrix organization
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(extra_cell_matrix)){
    if (Reactome_pathways1[i,2]==extra_cell_matrix[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Extracellular matrix organization"
    }
  }
}

#Loop para Gene expression (Transcription)
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(gene_expression)){
    if (Reactome_pathways1[i,2]==gene_expression[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Gene expression (Transcription)"
    }
  }
}

#Loop para Hemostasis
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(hemostasis)){
    if (Reactome_pathways1[i,2]==hemostasis[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Hemostasis"
    }
  }
}

#Loop para Immune System
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(immune_system)){
    if (Reactome_pathways1[i,2]==immune_system[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Immune System"
    }
  }
}

#Loop para Metabolism of proteins
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(metabolism_proteins)){
    if (Reactome_pathways1[i,2]==metabolism_proteins[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Metabolism of proteins"
    }
  }
}

#Loop para Metabolism of RNA
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(metabolism_rna)){
    if (Reactome_pathways1[i,2]==metabolism_rna[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Metabolism of RNA"
    }
  }
}

#Loop para Mitophagy
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(mitophagy)){
    if (Reactome_pathways1[i,2]==mitophagy[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Mitophagy"
    }
  }
}

#Loop para Muscle contraction
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(muscle_cons)){
    if (Reactome_pathways1[i,2]==muscle_cons[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Muscle contraction"
    }
  }
}

#Loop para Neuronal System
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(neuronal_system)){
    if (Reactome_pathways1[i,2]==neuronal_system[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Neuronal System"
    }
  }
}

#Loop para Organelle biogenesis and maintenance
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(organelle_biogenesis)){
    if (Reactome_pathways1[i,2]==organelle_biogenesis[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Organelle biogenesis and maintenance"
    }
  }
}

#Loop para DNA Repair
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(repairdna)){
    if (Reactome_pathways1[i,2]==repairdna[j,1]){
      Reactome_pathways1$Clasificador[i]<-"DNA Repair"
    }
  }
}

#Loop para DNA Replication
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(replication)){
    if (Reactome_pathways1[i,2]==replication[j,1]){
      Reactome_pathways1$Clasificador[i]<-"DNA Replication"
    }
  }
}

#Loop para Reproduction
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(reproduction)){
    if (Reactome_pathways1[i,2]==reproduction[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Reproduction"
    }
  }
}

#Loop para Transport of small molecules
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(transport_small_molecules)){
    if (Reactome_pathways1[i,2]==transport_small_molecules[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Transport of small molecules"
    }
  }
}

#Loop para Vesicle-mediated transport
for(i in 1:nrow(Reactome_pathways1)){
  for (j in 1:nrow(vesicle_transport)){
    if (Reactome_pathways1[i,2]==vesicle_transport[j,1]){
      Reactome_pathways1$Clasificador[i]<-"Vesicle-mediated transport"
    }
  }
}
