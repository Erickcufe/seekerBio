## code to prepare `populations_1000G` dataset goes here
populations_1000G <- completos[[1]]$population


usethis::use_data(populations_1000G, overwrite = TRUE)
