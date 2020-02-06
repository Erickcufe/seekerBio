## code to prepare `HouseKeeping` dataset goes here

HouseKeeping <- read.delim("data-raw/HouseKeeping_final.txt", header = FALSE)
colnames(HouseKeeping) <- c("Gen", "ID")

# usethis::use_data(HouseKeeping)
