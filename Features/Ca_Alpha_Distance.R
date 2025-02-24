# Implementation of the Alpha Carbon Distance feature

# Required libraries
library(bio3d)
library(magrittr)
library(dplyr)
library(tidyverse)

# Extract alpha carbons and calculate the Euclidean distance between wild-type and pathogenic pairs ----
# Input: native/mutated pdb file
get.ca.dist.mtx <- function(pdb) {
  
  ca_indices <- which(pdb$atom$elety == "CA")
  
  pdb$xyz <- matrix(pdb$xyz, ncol = 3, byrow = TRUE)
  
  coords <- pdb$xyz[ca_indices, , drop = FALSE]
  
  distance_matrix <- as.matrix(dist(coords))
  
  return(distance_matrix)
}

# Calculate the distance RMSD ----
# Input: path to the native and mutated pdb files
# Output: rmsd value of the difference between native and mutant mtx
dist.rmsd <- function(wt_pdb_path, mt_pdb_path) {
  
  wt_pdb <- read.pdb(wt_pdb_path)
  mt_pdb <- read.pdb(mt_pdb_path)
  
  # same pdb, will be zero
  dist_mtx_1 <- get.ca.dist.mtx(wt_pdb)
  dist_mtx_2 <- get.ca.dist.mtx(mt_pdb)
  
  if (!all(dim(dist_mtx_1) == dim(dist_mtx_2))) {
    stop("Structures have different number of C-alpha atoms!")
  }
  
  diff_mtx <- abs(dist_mtx_1 - dist_mtx_2)
  
  # rmsd function from bio3d package
  # rmsd <- rmsd(diff_mtx)
  rmsd <- sqrt(mean(diff_mtx^2))
  return(rmsd)
}
