# Calculates Delta SASA using the output from Dr.SASA

# Required Libraries
library(bio3d)
library(dplyr)
library(stringr)

# Normalize SASA function
calculate_sasa <- function(path_asa_file, pos) {
  
  # Load SASA file
  pdb_sasa <- read.csv(path_asa_file, skip = 1, header = FALSE, sep = "")

  # DR.SASA does not assign column names to the csv output - essential to include them
  colnames(pdb_sasa) <- c("ATOM", "Atom_Number", "Atom_Type", "Resid", "Chain", "Resno", "x", "y", "z", "occupancy", "SASA", "Type", "Element", "Random")
  
  # Filter for target residue position
  sasa_target_pos <- pdb_sasa %>% dplyr::filter(Resno == pos)
  
  # Calculate sum of SASA values
  sasa_target_sum <- sum(as.numeric(sasa_target_pos$SASA), na.rm = TRUE)
  
  return(sasa_target_sum)
}

# Calculate Delta SASA function
calculate_delta_sasa <- function(path_wt_pdb, path_mt_pdb, path_wt_asa, path_mt_asa, pos, max_sasa_value) {
  
  # Calculate SASA for wild_type and mutant structures
  wt_sasa <- calculate_sasa(path_wt_asa, pos)
  mt_sasa <- calculate_sasa(path_mt_asa, pos)
  
  # Handle cases where SASA is missing
  if (is.na(wt_sasa)) {
    return(list(wild_type_sasa = NA_real_, mutant_sasa = NA_real_, delta_sasa = NA_real_, delta_sasa_norm = NA_real_))
  }
  
  # Check if both SASA values are zero
  if (wt_sasa == 0 && mt_sasa == 0) {
    return(list(wild_type_sasa = 0, mutant_sasa = 0, delta_sasa = 0, delta_sasa_norm = 0))
  }
  
  # Compute delta SASA
  delta_sasa <- wt_sasa - mt_sasa
  
  # Normalize delta SASA using max SASA for the wild-type residue
  delta_sasa_norm <- ifelse(!is.na(max_sasa_value) && max_sasa_value > 0, delta_sasa / max_sasa_value, NA_real_)
  
  return(list(wild_type_sasa = wt_sasa, mutant_sasa = mt_sasa, delta_sasa = delta_sasa, delta_sasa_norm = delta_sasa_norm))
}

# Initialize results dataframe
results_df <- data.frame(
  uniprotid = character(),
  mutation = character(),
  filename = character(),
  position = numeric(),
  delta_sasa = numeric(),
  delta_sasa_norm = numeric(),
  stringsAsFactors = FALSE
)
