# Code to reproduce the results presented in my PhD Dissertation
This repository contains the necessary code to reproduce the results reported in my PhD dissertation tittled:

It includes scripts for:
- Predicting wild-type and mutated protein structures using ParaFold (HPC environment only)
- Source code for the developed structural features
- Source code to run the SIESTA-suite models: SIESTA, SIESTA-Str, and SIESTA-Seq

---

# Protein Structure Feature Extraction

This repository provides R scripts for extracting structural features from **wild-type and mutant** protein structures. These features help analyze the impact of mutations on protein stability and function.

## 📌 Features

Each feature is implemented as a standalone **R script**. The scripts can be found in the `Features` directory.


### **1. Alpha Carbon Distance**
- **Script:** `features/CA_Alpha_Distance.R`
- **Description:** Computes the distance between C-alpha atoms of wild-type and mutant structures.
- **Input:**
  - Wild-type PDB file  
  - Mutant PDB file  
- **Usage:**
  ```r
  dist.rmsd("path/to/wt.pdb", "path/to/mt.pdb")
  ```

### **2. Delta Alpha Carbon pLDDT**
- **Script:** `Features/Delta_CA_pLDDT.R`
- **Description:** Computes the difference between pLDDT values of C-alpha atoms of wild-type and mutant structures.
- **Input:**
  - Wild-type PDB file  
  - Mutant PDB file  
- **Usage:**
```r
   delta.ca.plddt.norm("path/to/wt.pdb", "path/to/mt.pdb")
  ```

### **3. Delta SASA**

- **Script:** `Features/Delta_SASA_Normalized_Max.R`
- **Description:** Computes the difference in SASA values.

- **Input:**  
  - Wild-type PDB file  
  - Mutant PDB file  
  - Wild-type ASA file  
  - Mutant ASA file  
  - Position where the variation occurs (int)  
  - Max SASA value (float)  

- **Usage**:
```r
calculate_delta_sasa("path/to/wt.pdb",
                     "path/to/mt.pdb",
                     "path/to/wt.asa",
                     "path/to/mt.asa",
                      52,
                      219,3),
```

### **4. MJ Potential Mutant**

- **Script:** `Features/MJ_Potential_Mutant.R`
- **Description:** Computes the MJ Potential Values between the identified contacts of the wild-type and mutated structures.
- **Input:**  
  - Wild-type avg_side_chain_table (csv)
  - Mutant avg_side_chain_table (csv)
  - MJ Table 
  - Wild-type amino acid residue (char)
  - Wild-type position (int)
  - Mutant amino acid residue (int) 
- **Usage**:
```r
calculate_dsum("path/to/wt_avg_side_chain_table.csv",
               "path/to/mt_avg_side_chain_table.csv",
               "path/to/mj_table_file.csv",
               "wt_residue",
                wt_location,
                mt_residue)
```


### **5. dRMS Local**

- **Script:** `Features/dRMS_Local.R`
- **Description:** Computes the MJ Potential Values between the identified contacts of the wild-type and mutated structures.
- **Input:**  
  - Wild-type JSON
  - Wild-type amino acid residue (char)
  - Mutant predicted structure in PDB format
  - Wild-type position (int)
- **Usage**:
```r
calculate_dsum("path/to/wt.json",
               "path/to/wt.pdb",
               "path/to/mt.pdb",
               "pos"
```





