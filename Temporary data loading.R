# Load needed packages
library(dplyr)
library(tidyr)
library(tibble)

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

# Extract columns (TFs) which have sufficiently large mean and variance








# Average over cell types
TF_averaged <- TF_array_filtered %>% 
  select(append(c("cell_ontology_class"), 
                tf_names)) %>% 
  group_by(cell_ontology_class) %>% 
  summarise(across(tf_names, mean, na.rm = TRUE)) %>% 
  column_to_rownames(var = "cell_ontology_class") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "tf") 

a <- TF_array_filtered %>% 
  mutate(across(where(is.numeric), ~ log1p(.))) %>% 
  














