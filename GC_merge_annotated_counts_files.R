# This script merges all annotated cell counts table into a summary cell counts table

# Load necessary package
library(dplyr)
library(tidyr)
library(readr)

# Read input parameters from the text file
input_file <- "Instructions_adj_GC_R2spleen2.txt"
params <- readLines(input_file)

# Parse parameters
for (param in params) {
  eval(parse(text = param))
}

# Set the working directory
setwd(working_directory)

# Get the list of all CSV files in the folder
files <- list.files(working_directory, pattern = "_counts\\.csv$", full.names = TRUE)

# Initialize an empty list to store the data and a set to track all unique cell types
data_list <- list()
all_cell_types <- character(0)  # Store all unique cell types encountered

# Loop over each file and load the data
for (file in files) {
  # Extract the unique part of the filename with regex (xxxx-xxxx format)
  filename <- basename(file)
  follicle_name <- sub("(^\\d{4}-\\d{4}).*", "\\1", filename)
  
  # Read the CSV file into a data frame
  df <- read.csv(file, header = FALSE, skip=1)
  
  # Track all unique cell types across files
  all_cell_types <- union(all_cell_types, df[[1]])
  
  # Convert the first column (cell category) into row names
  rownames(df) <- df[[1]]
  
  # Extract the counts column and convert to a named vector
  counts <- setNames(df[[2]], rownames(df))
  
  # Convert to a data frame with counts as one row
  counts_df <- as.data.frame(t(counts))
  
  # Add the follicle name as a column
  counts_df$Follicle <- follicle_name
  
  # Store the data frame in the list
  data_list[[follicle_name]] <- counts_df
}

# Combine all data frames, filling in missing values with 0
merged_annotated_counts <- bind_rows(data_list) %>%
  replace(is.na(.), 0)  # Replace any NAs with 0

# Move the 'Follicle' column to the first position
merged_annotated_counts <- merged_annotated_counts %>%
  select(Follicle, everything())


# Save the merged counts to a CSV file
write_csv(merged_annotated_counts, "merged_annotated_counts.csv")

