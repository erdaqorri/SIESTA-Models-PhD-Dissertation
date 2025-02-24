# Calculates the dRMS Local feature using a 13A distance filtering, it does not apply the pLDDT filtering

# Required libraries
library(jsonlite)
library(bio3d)
library(dplyr)
library(magrittr)
library(reshape2)
library(stringr)


process_protein_structures <- function(json_path, pdb_wild_type_path, mutated_pdb_path, position) {
  
  # Extract identified contacts from Arpeggio ----
  if (!file.exists(json_path)) {
    stop("JSON file does not exist: ", json_path)
  }
  if (!is.numeric(position) || position <= 0) {
    stop("Invalid position: ", position)
  }
  
  # Flatten the json file into a dataframe
  json_df <- jsonlite::fromJSON(json_path) %>%
    jsonlite::flatten() %>%
    as.data.frame()
  
  # Extract contacts where the wild_type residue is in the bgn pos
  contacts_bgn <- json_df %>% filter(bgn.auth_seq_id == position)
  # Extract contacts where the wild_type residue is in the end pos
  contacts_end <- json_df %>% filter(end.auth_seq_id == position)
  
  # Reverse the col names of the contacts made by the end res for easier joining of the contacts
  contacts_end <- contacts_end %>%
    rename(
      bgn.label_comp_id = `end.label_comp_id`,
      bgn.auth_seq_id = `end.auth_seq_id`,
      bgn.auth_atom_id = `end.auth_atom_id`,
      end.label_comp_id = `bgn.label_comp_id`,
      end.auth_seq_id = `bgn.auth_seq_id`,
      end.auth_atom_id = `bgn.auth_atom_id`
    )
  
  # Merge contacts identified in two directions together
  merged_contacts <- rbind(contacts_bgn, contacts_end)
  
  # Remove any duplicate entries
  merged_contacts_u <- unique(merged_contacts)
  
  # Rename the cols for easier processing
  merged_contacts_u <- merged_contacts_u %>%
    rename(
      res1 = `bgn.label_comp_id`,
      pos1 = `bgn.auth_seq_id`,
      atom_id1 = `bgn.auth_atom_id`,
      res2 = `end.label_comp_id`,
      pos2 = `end.auth_seq_id`,
      atom_id2 = `end.auth_atom_id`
    )
  
  # Interested only in the res, pos, and atom id of the interacting residues of the wild_type
  # Renamed to match the columns for extraction in the pdb file
  merged_contacts_u_f <- merged_contacts_u %>%
    rename(
      resid = res2,
      resno = pos2,
      elety = atom_id2
    )
  
  # Final list of contacts made by the wild_type residue in the mutated pos
  contacts_intial <- merged_contacts_u_f
  
  ###### Process the wild_type PDB file ######   
  wild_type_structure <- read.pdb(pdb_wild_type_path)
  
  # Extract only the heavy atoms, exclude the hydrogen atoms  
  heavy_atoms <- atom.select(wild_type_structure, "noh")
  coords <- wild_type_structure$atom[heavy_atoms$atom, ]
  # coords_filtered <- coords %>% select(elety, resid, resno, x, y, z)
  
  # Extract the xyz coords of the identified contacts based on the resid, resno, and elety
  filtered_pdb <- coords %>%
    inner_join(contacts_intial, by = c("resid", "resno", "elety"))
  
  # Adjacent pos do not give any additional information in term of contacts for the wild_type, and are removed
  # With two filters, -1, +1, and self contacts
  filtered_pdb_xyz <- filtered_pdb %>% select("resid", "resno", "elety", x, y, z)
  filtered_pdb_xyz <- filtered_pdb_xyz %>% filter(resno != (position - 1) &
                                                    resno != (position + 1) &
                                                    resno != position) %>%
    arrange(resno, resid, elety)
  
  # Create a wild_type distance matrix on the xyz coords of the identified contacts
  wt_dist_mtx <- as.matrix(dist(filtered_pdb_xyz[, c("x", "y", "z")]))
  
  # The rownames and colnames are lost from the matrix so add a new col to rename the rownames and colnames
  filtered_pdb_xyz$id_col <- paste(filtered_pdb_xyz$resid, filtered_pdb_xyz$elety, filtered_pdb_xyz$resno, sep = "_")
  
  rownames(wt_dist_mtx) <- filtered_pdb_xyz$id_col
  colnames(wt_dist_mtx) <- filtered_pdb_xyz$id_col
  
  # Check if the matrix is symmetric
  # If yes, extract the upper part of the matrix
  is_symmetric <- all.equal(wt_dist_mtx, t(wt_dist_mtx))
  if (is_symmetric != TRUE) {
    stop("The distance matrix is not symmetric. Exiting.")
  }
  
  # Extract the upper part of the matrix
  # Set NA the lower part of the matrix
  upp_t <- wt_dist_mtx
  upp_t[!upper.tri(wt_dist_mtx)] <- NA
  
  # Matrix to df for easier processing  
  wt_dist_mtx_df <- melt(upp_t, varnames = c("mtx_rows", "mtx_cols"), value.name = "dist")
  
  # Remove the NA and atomic contacts within a distance greater than 5
  wt_dist_mtx_df_filtered <- wt_dist_mtx_df %>% filter(!is.na(dist),
                                                       dist !=0,
                                                       dist <= 13)
  
  # Gly tends to form a small number of contacts so leads to zero contacts after the dist filtering
  # Set drms_local to 0 if Gly does not have any contacts
  if (nrow(wt_dist_mtx_df_filtered) == 0) {
    print("No filtered contacts found. Setting drms_local to 0.")
    drms_local <- 0
    
    # Return early with a modified result structure
    return(list(
      wild_type_results = list(contacts_initial = contacts_intial, contacts_final = data.frame(), filtered_uniq_df = data.frame()),
      mutated_results = list(
        mutated_contacts = data.frame(),
        mutated_distance_matrix = matrix(0),
        matching_rows = data.frame(),
        drms_int = drms_local)
    ))
  }
  
  # Identify unique contacts and remove the reversed ones
  wt_dist_mtx_df_filtered_uniq <- unique(wt_dist_mtx_df_filtered)
  
  # Every identified contact might make contacts with different atoms
  # Extract the unique ones from the rows and cols and merge them into a vec and create a df
  wt_dist_mtx_df_filtered_uniq$mtx_rows %>% unique() -> rows_uniq
  wt_dist_mtx_df_filtered_uniq$mtx_cols %>% unique() -> cols_uniq
  
  rows_uniq <- as.vector(rows_uniq)
  cols_uniq <- as.vector(cols_uniq)
  merged_vec <- c(rows_uniq, cols_uniq)
  merged_vec_uniq <- unique(merged_vec)
  df <- data.frame(unique_contacts = merged_vec_uniq)
  
  # Extract the list of identified contacts within <= 5 to search in the mutant structure
  contacts_df <- df %>%
    mutate(resid = sub("([A-Z]{3})_.*", "\\1", unique_contacts),
           resno = sub(".*_(\\d+)$", "\\1", unique_contacts),
           elety = sub("^[A-Z]{3}_(.*)_[0-9]+$", "\\1", unique_contacts))
  
  # Col resno by default gets set as char and raises errors so change it into int
  contacts_df$resno <- as.integer(contacts_df$resno)
  
  ###### Process the mutated PDB file ######   
  mutated_structure <- read.pdb(mutated_pdb_path)
  heavy_atoms_mt <- atom.select(mutated_structure, "noh")
  mt_pdb_xyz <- mutated_structure$atom[heavy_atoms_mt$atom, ]
  # mt_pdb_xyz <- coords_mt %>% select("resid", "resno", "elety", x, y, z)
  
  # Extract the xyz of the contacts identified in the wild_type structure
  mt_contacts_xyz <- mt_pdb_xyz %>%
    inner_join(contacts_df, by = c("resid", "resno", "elety"))
  
  mt_contacts_xyz <- mt_contacts_xyz %>%
    arrange(resno, resid, elety)
  
  # Add a colname to matrix to identify the interacting atom pairs
  mt_contacts_xyz$combined <- paste(mt_contacts_xyz$resid, mt_contacts_xyz$elety, mt_contacts_xyz$resno, sep = "_")
  
  # Create the mutant distance matrix on these xyz of the identified contacts in the mutant
  # To check for changes in their geometry
  mt_dist_mtx <- as.matrix(dist(mt_contacts_xyz[, c("x", "y", "z")]))
  
  rownames(mt_dist_mtx) <- mt_contacts_xyz$combined
  colnames(mt_dist_mtx) <- mt_contacts_xyz$combined
  
  is_symmetric <- all.equal(mt_dist_mtx, t(mt_dist_mtx))
  if (is_symmetric != TRUE) {
    stop("The distance matrix is not symmetric. Exiting.")
  }
  
  upp_t_mt <- mt_dist_mtx
  upp_t_mt[!upper.tri(mt_dist_mtx)] <- NA
  
  mt_dist_mtx_tb <- melt(upp_t_mt, varnames = c("mtx_rows_mt", "mtx_cols_mt"), value.name = "dist")
  mt_dist_mtx_tb_filtered <- mt_dist_mtx_tb %>% filter(!is.na(dist),
                                                       dist != 0)
  
  mt_dist_mtx_tb_filtered_uniq <- unique(mt_dist_mtx_tb_filtered)
  
  # Extract those pairs from the mutant unique interactions that are present in the wild_type one 
  mt_matching_rows <- semi_join(mt_dist_mtx_tb_filtered_uniq, wt_dist_mtx_df_filtered_uniq, by = c("mtx_rows_mt" = "mtx_rows", "mtx_cols_mt" = "mtx_cols"))
  
  # Rename the dist into dist_mt
  wt_dist_mtx_contacts <- wt_dist_mtx_df_filtered_uniq %>% rename("dist_wt" = dist)
  mt_matching_rows <- mt_matching_rows %>% rename("dist_mt" = dist)
  
  wt_mt_results <- wt_dist_mtx_contacts %>%
    left_join(mt_matching_rows, by = c("mtx_rows" = "mtx_rows_mt", "mtx_cols" = "mtx_cols_mt")) %>% 
    mutate(`wt-mt` = (dist_wt - dist_mt)^2)
  
  # apply the formula
  drms_local <- sqrt(sum(wt_mt_results$`wt-mt`)) / nrow(wt_dist_mtx_contacts)
  
  return(list(
    wild_type_results = list(contacts_initial = contacts_intial, contacts_final = contacts_df, filtered_uniq_df = wt_dist_mtx_df_filtered_uniq),
    mutated_results = list(
      mutated_contacts = mt_contacts_xyz,
      mutated_distance_matrix = mt_dist_mtx,
      matching_rows = mt_matching_rows,
      drms_int = drms_local)
  ))
}
