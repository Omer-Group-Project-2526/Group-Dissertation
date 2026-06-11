# Possible to-dos
# 1. Record cell numbers of each type
# 2. Record mean and std for each TF (so a better assessment can be made on cutoff)
# 3. Considerations on exact method of normalisation (expect this to have a small impact)

# Load needed packages
library(dplyr)
library(tidyr)
library(tibble)

# Set working directory (can be done manually or with code)
if (getwd() != `DESIRED PATH`) {
  setwd(`DESIRED PATH`)
}

normalisation_method <- menu(c("L1 (Manhattan)", "L2 (Euclidean)"), title="Choose normalisation methodology:")

# Load annotation file needed for the translation from codes to cell types
FACS_annotation <- read.csv2("00_facs_raw_data/annotations_FACS.csv", sep = ",") %>% 
  as.data.frame() %>% 
  extract(cell, into = c("cell.id", "plate.barcode", "mouse.id"),
          regex = "^([^.]+)\\.([^.]+)\\.([^.]+)\\.",
          remove = FALSE) %>% 
  filter(cell_ontology_class != "unknown")

# List of transcription factors (filter down from proteins)
TF_list <- read.csv2("Mus_musculus_TF.txt", sep = "\t", strip.white = TRUE)$Symbol

# Find all FACS files to load cell and TF data from (stored in separate csvs)
TF_csvs <- list.files("00_facs_raw_data/FACS")

# Initialise array
TF_array <- data.frame()

# Iterate over files, filter for the TFs and then append to larger df
for (indx in 1:length(TF_csvs)) {
    print(paste(indx, TF_csvs[indx], sep = ": "))
  
    TF_to_append <- as.data.frame(read.csv2(paste("00_facs_raw_data/FACS", TF_csvs[indx], sep = "/"), sep = ",")) %>% 
      rename(tf = X) %>% 
      filter(tf %in% TF_list)
    
    if(indx == 1) {
      TF_array <- TF_to_append
    } else {
      TF_array <- full_join(TF_array, TF_to_append, by = "tf")
    }
    
    print(" -> Appended")
}

# Reorder for easier usage
TF_array_filtered <- TF_array %>% 
  column_to_rownames(var = "tf") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "cell") %>% 
  inner_join(FACS_annotation, by = join_by(cell == cell)) %>% 
  mutate(across(where(is.numeric), ~ log1p(.)))

# Non-TF column names
non_tf_names <- setdiff(colnames(TF_array_filtered),
                        TF_list)

# TF column names
tf_names <- intersect(colnames(TF_array_filtered),
                      TF_list)

# Now filter TFs to have sufficient mean and variance
TF_array_filtered <- bind_cols(
  TF_array_filtered %>% 
    select(all_of(non_tf_names)),
  TF_array_filtered %>% 
    select(all_of(tf_names)) %>% 
    select(where(~ mean(., na.rm = TRUE) > log(4) & 
                   sd(., na.rm = TRUE)   > log(4)))
)

# Update TF names as now many have been filtered/removed
tf_names_filtered <- intersect(colnames(TF_array_filtered),
                               tf_names)

# Average over cell types
TF_averaged <- TF_array_filtered %>% 
  select(append(c("cell_ontology_class"), 
                tf_names_filtered)) %>% 
  group_by(cell_ontology_class) %>% 
  summarise(across(tf_names_filtered, mean, na.rm = TRUE)) %>% 
  column_to_rownames(var = "cell_ontology_class") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "tf")

# Normalise averaged data - this is Xi (or at least the internal data is!)
TF_normalised <- TF_averaged %>% 
  column_to_rownames(var = "tf") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "cell_ontology_class") %>% 
  mutate(across(all_of(tf_names_filtered), 
                if (normalisation_method == 1) {
                  ~ . / sum(., na.rm = TRUE)
                } else {
                  ~ . / sqrt(sum(. ^ 2, na.rm = TRUE))
                }))

# Write Xi to a CSV to be saved and potentially shared
write.csv(TF_normalised, "FACS_Xi_matrix.csv", sep = ",")
