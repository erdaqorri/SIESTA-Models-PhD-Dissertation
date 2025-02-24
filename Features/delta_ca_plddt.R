# Calculate the difference in pLDDT scores assigned to each carbon alpha between wild type and mutated protein structure
# Does not apply per-residue pLDDT filtering (a global one was applied on the structures)

library(bio3d)
library(dplyr)

delta.ca.plddt.norm <- function(wt_pdb, mt_pdb){
  
  # Read in the wild_type and mutant pdb files
  wild_type_pdb <- read.pdb(wt_pdb)
  mutant_pdb <- read.pdb(mt_pdb)
  
  # Extract the CA carbons from the wild_type structure
  wild_type.atom.df <- wild_type_pdb$atom
  wild_type.atom.df <- wild_type.atom.df %>% mutate(elety = ifelse(elety == "CA",
                                                             paste0("CA_", cumsum(elety == "CA")),
                                                             elety))
  
  # Extract only the CA atoms
  wild_type.ca <- wild_type.atom.df[grepl("CA_", wild_type.atom.df$elety), ]
  
  if (nrow(wild_type.ca) == 0) {
    stop("Attention! No CA atoms where found in the input PDB file! Check it and return!")
  }
  
  # Extract the CA atoms from the mutated structure
  mutant.atom.df <- mutant_pdb$atom
  mutant.atom.df <- mutant.atom.df %>% mutate(elety = ifelse(elety == "CA",
                                                             paste0("CA_", cumsum(elety == "CA")),
                                                             elety))
  
  mutant.ca <- mutant.atom.df[grepl("CA_", mutant.atom.df$elety), ]
  
  if (!length(wild_type.ca$elety) == length(mutant.ca$elety)) {
    stop("The number of CA atoms is not the same!")
  }
  
  if (!identical(wild_type.ca$elety, mutant.ca$elety)) {
    stop("The CA information is corrupt!")
  }
  
  # Merge the wild_type and mutant dataframes
  # Interested only in the elety, resid, and b columns
  wild_type.ca_subset <- wild_type.ca %>% select(elety, resid, b)
  wild_type.ca_subset$order <- 1:nrow(wild_type.ca_subset)
  
  wt.mt.merged <- merge(wild_type.ca_subset, mutant.ca[, c("elety", "b", "resid")],
                        by = "elety",
                        suffixes = c("_wild_type", "_mutant"))
  
  wt.mt.merged_ord <- wt.mt.merged[order(wt.mt.merged$order), ]
  
  # Remove the temporary 'order' column
  wt.mt.merged_ord$order <- NULL
  
  # substract them
  wt.mt.pldt <- wt.mt.merged_ord %>% mutate(wt_mt_diff = abs(b_wild_type - b_mutant))
  wt.mt.pldt.sum <- wt.mt.pldt$wt_mt_diff %>% sum()
  final_sum <- wt.mt.pldt.sum / length(wt.mt.pldt$elety)
  
  return(final_sum)
}
