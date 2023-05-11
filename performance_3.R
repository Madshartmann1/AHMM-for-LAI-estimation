# Get file name from command line argument
input_files <- commandArgs(trailingOnly = TRUE)

# Load files
mutation_list <- data.frame(read.table("mutations.list"))
theoretical <- data.frame(read.table("theoretical_intervals.tsv", header = TRUE))



for (input_file in input_files){
  
  AHMM_est <- data.frame(read.table(input_file, header = TRUE))
  
  #setting up counters and exeptions
  Match <- Partial <- Mismatch <- ancient_mut <- 0
  ancient_puplations <- c("ancestral", "AMH", "OOA", "ADMIX")
  no_punishment <- c("CHB","KAR")
  print(input_file)
  #Main comparison code
  for (p in mutation_list[,1]) {
#    if (p %% 100 == 0) {
#      print(paste("Analyzing position:", p))
#    }
    if (p < AHMM_est[,2]){
      next }
    #if position is at start of interval
    if (p %in% AHMM_est$start) {
      in_theoretical <- p >= theoretical$start & p <= theoretical$end
      if (any(in_theoretical)) {
        theoretical_name <- theoretical$Ancestral_population[in_theoretical]
        AHMM_name <- AHMM_est$mutation_ancestery_names[p == AHMM_est$start]
#        if (is.na(theoretical_name) || is.na(AHMM_name)) {
#          print("Names missing")
#          break
#        } 
        if (!(theoretical_name %in% ancient_puplations)) {
          if (theoretical_name == AHMM_name) {
            Match <- Match + 1
          } else if (theoretical_name %in% no_punishment & AHMM_name %in% no_punishment) {
            Match <- Match + 1
            
          } else if (grepl("/", AHMM_name) == T) {
            if (theoretical_name == sapply(strsplit(AHMM_name, "/"), `[`, 1) || theoretical_name == sapply(strsplit(AHMM_name, "/"), `[`, 2)) {
              Partial <- Partial + 1
            } else if (theoretical_name %in% no_punishment & sapply(strsplit(AHMM_name, "/"), `[`, 1) %in% no_punishment |
                       theoretical_name %in% no_punishment & sapply(strsplit(AHMM_name, "/"), `[`, 2) %in% no_punishment) {
              Partial <- Partial + 1
              
            } else {
              Mismatch <- Mismatch + 1
            }
          } else {
            Mismatch <- Mismatch + 1
          }
        } 
        else if  (theoretical_name %in% ancient_puplations) {
          ancient_mut <- ancient_mut + 1} 
      } else {
        Mismatch <- Mismatch + 1
      }
    } 
    #if position is at end of interval
    else if (p %in% AHMM_est$end) {
      in_theoretical <- p >= theoretical$start & p <= theoretical$end
      if (any(in_theoretical)) {
        theoretical_name <- theoretical$Ancestral_population[in_theoretical]
        AHMM_name <- AHMM_est$mutation_ancestery_names[p == AHMM_est$end]
#        if (is.na(theoretical_name) || is.na(AHMM_name)) {
#          print("Names missing")
#          break
#        } 
        if (!(theoretical_name %in% ancient_puplations)) {
          if (theoretical_name == AHMM_name) {
            Match <- Match + 1
          }
          else if (theoretical_name %in% no_punishment & AHMM_name %in% no_punishment) {
            Match <- Match + 1
            
          }
          else if (grepl("/", AHMM_name) == T) {
            if (theoretical_name == sapply(strsplit(AHMM_name, "/"), `[`, 1) || theoretical_name == sapply(strsplit(AHMM_name, "/"), `[`, 2)) {
              Partial <- Partial + 1
            } else if (theoretical_name %in% no_punishment & sapply(strsplit(AHMM_name, "/"), `[`, 1) %in% no_punishment |
                       theoretical_name %in% no_punishment & sapply(strsplit(AHMM_name, "/"), `[`, 2) %in% no_punishment){
              Partial <- Partial + 1
            }
            else {
              Mismatch <- Mismatch + 1
            }
          }
          else {
            Mismatch <- Mismatch + 1}
        } 
        else if  (theoretical_name %in% ancient_puplations) {
          ancient_mut <- ancient_mut + 1} 
        else {
          Mismatch <- Mismatch + 1
        }
      }
    } 
    #if position is in interval
    # else if (p < min(theoretical$start) | p > max(theoretical$end)) {
    #   print("Position not found in theoretical intervals")
    #   break
    # } 
    # else if (p < min(AHMM_est$start) | p > max(AHMM_est$end)) {
    #   print("Position not found in AHMM_est intervals")
    #   break
    # } 
    else if (!(p < min(theoretical$start) | p > max(theoretical$end) & p < min(AHMM_est$start) | p > max(AHMM_est$end))){
      j <- which(sapply(1:nrow(theoretical), function(j) p >= theoretical$start[j] & p <= theoretical$end[j]))
      k <- which(sapply(1:nrow(AHMM_est), function(k) p >= AHMM_est$start[k] & p <= AHMM_est$end[k]))
      theoretical_name <- theoretical$Ancestral_population[j]
      AHMM_name <- AHMM_est$mutation_ancestery_names[k]
#      if (is.na(theoretical_name) | is.na(AHMM_name)) {
#        print("Interval Names missing")
#        break }
       if (theoretical_name %in% ancient_puplations) {
        ancient_mut <- ancient_mut + 1 }
       else {
        if (theoretical_name == AHMM_name) {
          Match <- Match + 1 }
        else if (theoretical_name %in% no_punishment & AHMM_name %in% no_punishment) {
          Match <- Match + 1 }
        else if (grepl("/", AHMM_name) == TRUE) {
          AHMM_names <- unlist(strsplit(AHMM_name, "/"))
          if (theoretical_name %in% AHMM_names | theoretical_name %in% no_punishment & any(AHMM_names %in% no_punishment)) {
            Partial <- Partial + 1 }
           else {
            Mismatch <- Mismatch + 1}
          }
         else {
          Mismatch <- Mismatch + 1
        }
      }
    }
#    else {
#      print("Fallthrough")}
  }                                                                  
  
  
  # Save results
  result_df <- data.frame("Match" = Match,"Partial" = Partial ,"Mismatch" = Mismatch, "Ancients" = ancient_mut, "Total" = Match + Mismatch + Partial + ancient_mut)
  print(result_df)
  
  # Get output filename and write CSV
  output_filename <- paste0("performance/", sub("_intervals.tsv", ".csv", basename(input_file)))
  write.csv(result_df, output_filename)
}

