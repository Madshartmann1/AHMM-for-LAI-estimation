library(dplyr)
library(tidyr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script_name.R input_file output_file")
}
input_file <- args[1]
output_file <- args[2]

theoretical_df = read.table(input_file, header = T)

theoretical_intervals <- theoretical_df %>%
  mutate(
    change_category = Ancestral_population != lag(Ancestral_population, default = first(Ancestral_population)),
    group = cumsum(change_category)
  ) %>%
  group_by(Ancestral_population, group) %>%
  summarize(
    start = first(Position),
    end = last(Position)
  ) %>%
  ungroup() %>%
  select(-group) %>%
  arrange(start)


#save output as table 
write.table(theoretical_intervals, file = output_file, sep = "\t", quote = F, row.names = F)
