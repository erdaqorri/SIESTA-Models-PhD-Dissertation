# MJ-Potential of the Mutant Script based on the MJ-Potential feature from Ozkan et al, 2024
# Requires preprocessed output of the cif files using the initial implementation described in Ozkan et al, 2024 (attached as a script)

# Required libraries
library(tidyverse)


calculate_dsum <- function(wild_type_file, mutated_file, mj_table_file, wt_residue, wt_location, mt_residue) {
  # Read wild_type and mutated contact data
  wild_type_contacts <- read.csv(wild_type_file)
  mutated_contacts <- read.csv(mutated_file)
  
  # Read MJ potential table
  mj_table <- read.csv(mj_table_file, row.names = 1, check.names = FALSE)
  
  # Define function to get MJ potential
  get_mj_potential <- function(res1, res2) {
    value <- mj_table[res1, res2]
    if (is.na(value)) {
      value <- mj_table[res2, res1]
    }
    return(value)
  }
  
  # Filter wild_type contacts based on WT residue and location
  filtered_wild_type <- wild_type_contacts %>%
    filter(wt == wt_residue, pos == wt_location) %>%
    select(residue, wt, pos, contact_residue, contact_wt)
  
  # Calculate MJ sum for wild_type contacts (default to 0 if empty)
  if (nrow(filtered_wild_type) > 0) {
    filtered_wild_type$mj_pair <- mapply(get_mj_potential, filtered_wild_type$wt, filtered_wild_type$contact_wt)
    wild_type_mj_sum <- sum(filtered_wild_type$mj_pair, na.rm = TRUE)
  } else {
    wild_type_mj_sum <- 0
  }
  
  # Filter mutated contacts based on MT residue and location
  filtered_mutated <- mutated_contacts %>%
    filter(wt == mt_residue, pos == wt_location) %>%
    select(residue, wt, pos, contact_residue, contact_wt)
  
  # Calculate MJ sum for mutated contacts (default to 0 if empty)
  if (nrow(filtered_mutated) > 0) {
    filtered_mutated$mj_pair_MT <- mapply(get_mj_potential, filtered_mutated$wt, filtered_mutated$contact_wt)
    mutant_mj_sum <- sum(filtered_mutated$mj_pair_MT, na.rm = TRUE)
  } else {
    mutant_mj_sum <- 0
  }
  
  # Calculate the difference in MJ potentials
  dsum <- mutant_mj_sum - wild_type_mj_sum
  
  return(dsum)
}
