library(dplyr)
library(tidyr)


#function for translating the IDs of AHMM to their real names
#only works if the simulation is made in the following order "YRI", "CEU", "CHB", "KAR"
translate_id <- function(id) {
  # extract the relevant part of the ID string
  id_part <- sub("^X", "", id)
  
  # split the ID part into its four components
  components <- strsplit(id_part, "\\.")[[1]]
  if (length(components) != 4) stop("Invalid ID format: length")
  
  x <- as.numeric(components[1])
  y <- as.numeric(components[2])
  z <- as.numeric(components[3])
  v <- as.numeric(components[4])
  
  names <- c("YRI", "CEU", "CHB", "KAR")
  
  if (x == 2) {
    return(names[1])
  } else if (x == 1) {
    if (y == 1) return(paste(names[1], names[2], sep = "/"))
    if (z == 1) return(paste(names[1], names[3], sep = "/"))
    if (v == 1) return(paste(names[1], names[4], sep = "/"))
  } else if (x == 0) {
    if (y == 2) return(names[2])
    if (z == 2) return(names[3])
    if (v == 2) return(names[4])
    if (y == 1 && z == 1) return(paste(names[2], names[3], sep = "/"))
    if (y == 1 && v == 1) return(paste(names[2], names[4], sep = "/"))
    if (z == 1 && v == 1) return(paste(names[3], names[4], sep = "/"))
  }
  
  stop("Invalid ID format: translate")
}




# Get input file names from command line arguments
input_files <- commandArgs(trailingOnly = TRUE)

# Loop through input files
for (input_file in input_files) {

  # Read the data
  df <- read.table(input_file, sep = "\t", header = T)
  print(input_file)
  #the chromosome number is not needed
  reduced <- na.omit(select(df,-chrom))
  #taking out the stuff we actually need 
  temp = select(reduced,-position)

  #making a max column that indicates which pop the position is most proberble to stem from 
  reduced$p_max_population <- colnames(temp)[apply(temp, 1, which.max)]
  #getting positions and populations
  position_df <- select(reduced,position,p_max_population)

  # calculate the ranges of positions for each population
  postion_intervals <- position_df %>%
    mutate(
      change_category = p_max_population != lag(p_max_population, default = first(p_max_population)),
      group = cumsum(change_category)
    ) %>%
    group_by(p_max_population, group) %>%
    summarize(
      start = first(position),
      end = last(position)
    ) %>%
    ungroup() %>%
    select(-group) %>%
    arrange(start)




  #translate the names to the real ones 
  translate_names <- postion_intervals$p_max_population
  mutation_ancestery_names <- sapply(translate_names,translate_id)
  #ready to be saved 
  mutation_ancestery_est <- data.frame(mutation_ancestery_names, start = postion_intervals$start, end = postion_intervals$end)

  # Construct output file name
  output_file <- paste0("interval_predictions/",sub(".posterior", "_intervals.tsv",basename(input_file)))

  # Save output to file
  write.table(mutation_ancestery_est, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}


  
