# Code to reproduce the results presented in my PhD Dissertation
This repository contains the necessary code to reproduce the results reported in my PhD dissertation tittled:

It includes scripts for:
- Predicting wild-type and mutated protein structures using ParaFold (HPC environment only)
- Source code for the developed structural features
- Source code to run the SIESTA model
---

# Protein Structure Feature Extraction

## Input Data
The required input data to test the features along with the feature-specific files can be found in the **Data/Input_PDB** directory.
This directory has the following structure:
1. Data/Input_PDB : contains three PDB files (wt, pathogenic, and benign)
     - P68871_wt_rlx_model.pdb
     - P68871_V127G_pathogenic_rlx_model.pdb
     - P68871_K133N_benign_rlx_model.pdb
2. Data/Delta_Sasa_Normalized_Max : contains three ASA files (wt, pathogenic, and benign) as well as Max SASA value that a residue can have
     - P68871_wt_rlx_model.asa
     - P68871_V127G_pathogenic_rlx_model.asa
     - P68871_K133N_benign_rlx_model.asa
     - max_value_sasa_per_amino_acid.csv

  The ASA files were generated using the *dr-sasa* algorithm described in *Judemir Ribeiro, Carlos Ríos-Vera, Francisco Melo, Andreas Schüller, Calculation of accurate interatomic contact surface areas for the quantitative analysis of non-bonded molecular interactions, Bioinformatics, Volume 35, Issue 18, September 2019, Pages 3499–3501, https://doi.org/10.1093/bioinformatics/btz062*. To generate the files provided here, dr-sasa was run on **mode 0** with default parameters.

The maximum SASA value of every amino acid residue was calculated on a curated collection of experimentally resolved protein structures, provided by de la Cruz lab.

3. Data/MJ_Potential_Mutant: contains three csv files (wt, pathogenic, benign) and the MJ energy potentials table
     - P68871_wt_contact_avg_side_chain_table.csv
     - P68871_V127G_pathogenic_contact_avg_side_chains.csv
     - P68871_K133N_benign_contact_avg_side_chains.csv
     - MJ_Potential_Table.csv
  
   The contact average side chains were calculated used a script provided by de la Cruz lab, while the MJ Potential Table was also provided by this lab.

4. Data/dRMS_Local: contains one JSON file of the wild type
     - P68871_wt_rlx_model.json
  
   The JSON file was generated using the Arpeggio algorithm (https://pdbeurope.github.io/api-webinars/webinars/web5/arpeggio.html)


## 📌 Features

Each feature is implemented as a standalone **R script**. The scripts can be found in the `Features` directory.


### **1. Alpha Carbon Distance**

Input PDB files can be retrieved from: 

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

## 📌 Run SIESTA

SIESTA is an XGBoost model trained on both sequence and structure derived features from the AlphaFold2 generated structures.

- **Script:** `SIESTA.py`
- **Input:**  
  - --features.csv: csv file containing the extracted features
  - --output_model: trained SIESTA model saved in pkl format
- **Usage**:

```bash
python3 SIESTA.py  --features_csv "siesta_feature_test.csv"  --output_model "siesta.pkl" 
```




