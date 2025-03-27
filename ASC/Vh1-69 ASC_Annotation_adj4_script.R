# Load necessary libraries
library(readr)
library(stringr)
library(dplyr)
library(dbscan)
library(FNN)
library(ggplot2)
library(imager)
library(png)

# Read input parameters
input_file <- "Instructions_adj_ASC_R1spleen1.txt"
params <- readLines(input_file)

# Parse parameters
for (param in params) {
  eval(parse(text = param))
}


# Set the working directory
setwd(working_directory)

# Get all relevant files in the working directory
files <- list.files(working_directory)

# Extract unique sample prefixes
sample_prefixes <- unique(str_extract(files, "^[0-9]{4}-[0-9]{4}"))

# Initialize lists to hold intensity data
labels <- c()
CD45_1_intensities <- c()
CD45_2_intensities <- c()

# Collect intensity data from all samples
for (prefix in sample_prefixes) {
  CD45_1_file <- paste0(prefix, "_CD45_1.csv")
  CD45_2_file <- paste0(prefix, "_CD45_2.csv")
  if (all(c(CD45_1_file, CD45_2_file) %in% files)) {
    CD45_1_props <- read_csv(CD45_1_file, show_col_types = FALSE)
    CD45_2_props <- read_csv(CD45_2_file, show_col_types = FALSE)
    sample_labels <- paste0(prefix, "_", CD45_1_props$label)
    labels <- c(labels, sample_labels)
    CD45_1_intensities <- c(CD45_1_intensities, CD45_1_props$mean_intensity)
    CD45_2_intensities <- c(CD45_2_intensities, CD45_2_props$mean_intensity)
    rm(CD45_1_props, CD45_2_props)
  }
}

# Combine data into a data frame
combined_intensities <- data.frame(
  label = labels,
  cd45_1_intensity = CD45_1_intensities,
  cd45_2_intensity = CD45_2_intensities
)

# Scale the data for DBSCAN
scaled_data <- combined_intensities %>%
  mutate(
    cd45_1_scaled = scale(cd45_1_intensity),
    cd45_2_scaled = scale(cd45_2_intensity)
  )

# Calculate dynamic `eps` based on k-nearest neighbor distances
k <- 5
knn_distances <- FNN::get.knn(scaled_data[, c("cd45_1_scaled", "cd45_2_scaled")], k = k)$nn.dist
dynamic_eps <- median(knn_distances[, k])

# Apply DBSCAN with dynamic `eps`
dbscan_result <- dbscan::dbscan(scaled_data[, c("cd45_1_scaled", "cd45_2_scaled")], eps = dynamic_eps, minPts = k)

# Append DBSCAN cluster assignments to data frame
combined_intensities$cluster <- dbscan_result$cluster

# Assign clusters to categories based on thresholds
thresholds <- list(
  WT = list(cd45_1_intensity = quantile(combined_intensities$cd45_1_intensity, 0.75), cd45_2_intensity = quantile(combined_intensities$cd45_2_intensity, 0.75)),
  VH169 = list(cd45_1_intensity = quantile(combined_intensities$cd45_1_intensity, 0.25), cd45_2_intensity = quantile(combined_intensities$cd45_2_intensity, 0.75)),
  "564" = list(cd45_1_intensity = quantile(combined_intensities$cd45_1_intensity, 0.75), cd45_2_intensity = quantile(combined_intensities$cd45_2_intensity, 0.25))
)
combined_intensities <- combined_intensities %>%
  mutate(classification = case_when(
    cluster == 0 ~ "Noise",
    abs(cd45_1_intensity - thresholds$WT$cd45_1_intensity) + abs(cd45_2_intensity - thresholds$WT$cd45_2_intensity) < dynamic_eps ~ "WT",
    abs(cd45_1_intensity - thresholds$VH169$cd45_1_intensity) + abs(cd45_2_intensity - thresholds$VH169$cd45_2_intensity) < dynamic_eps ~ "VH169",
    abs(cd45_1_intensity - thresholds$`564`$cd45_1_intensity) + abs(cd45_2_intensity - thresholds$`564`$cd45_2_intensity) < dynamic_eps ~ "564",
    TRUE ~ "Unclassified"
  ))

# Ensure no unclassified cells or noise
unclassified_indices <- which(combined_intensities$classification %in% c("Unclassified", "Noise"))
for (i in unclassified_indices) {
  cell <- combined_intensities[i, ]
  
  # Calculate distances to all thresholds
  distances <- sapply(names(thresholds), function(category) {
    abs(cell$cd45_1_intensity - thresholds[[category]]$cd45_1_intensity) +
      abs(cell$cd45_2_intensity - thresholds[[category]]$cd45_2_intensity)
  })
  
  # Reassign to the closest cluster
  closest_category <- names(which.min(distances))
  combined_intensities$classification[i] <- closest_category
}

# Remove suffixes to ensure consistency
combined_intensities <- combined_intensities %>%
  mutate(classification = gsub("\\.\\d+%?$", "", classification))

# Define a tolerance for desired percentages
tolerance <- 0.5

# Calculate desired counts based on percentages with tolerance
total_cells <- nrow(combined_intensities)
desired_counts <- round(total_cells * c(WT = percentage_B6, VH169 = percentage_Vh169, "564" = percentage_564))
allowed_ranges <- list()
for (category in names(desired_counts)) {
  lower_bound <- max(0, desired_counts[category] - round(tolerance * desired_counts[category]))
  upper_bound <- desired_counts[category] + round(tolerance * desired_counts[category])
  allowed_ranges[[category]] <- c(lower_bound, upper_bound)
}

# Redistribution function with tolerance
redistribute_with_tolerance <- function(data, thresholds, desired_counts, allowed_ranges) {
  initial_counts <- table(factor(data$classification, levels = names(desired_counts)))
  deviation <- desired_counts - initial_counts
  for (category in names(deviation)) {
    if (initial_counts[category] < allowed_ranges[[category]][1]) {
      deficit <- allowed_ranges[[category]][1] - initial_counts[category]
      
      # Reassign from the closest excess category
      excess_category <- names(which.max(deviation[deviation > 0]))
      if (!is.na(excess_category)) {
        excess_cells <- data %>%
          filter(classification == excess_category) %>%
          mutate(
            distance = abs(cd45_1_intensity - thresholds[[category]]$cd45_1_intensity) +
                       abs(cd45_2_intensity - thresholds[[category]]$cd45_2_intensity)
          ) %>%
          arrange(distance)
        
        # Assign deficit number of cells
        reassigned_cells <- head(excess_cells, deficit)
        data$classification[data$label %in% reassigned_cells$label] <- category
        deviation[category] <- deviation[category] + deficit
        deviation[excess_category] <- deviation[excess_category] - deficit
      }
    }
  }
  return(data)
}

# Apply redistribution logic
combined_intensities <- redistribute_with_tolerance(
  combined_intensities,
  thresholds,
  desired_counts,
  allowed_ranges
)

# Second redistribution
min_564 <- round(desired_counts["564"] * 0.85)
min_WT <- round(desired_counts["WT"] * 1.5)
min_VH169 <- round(desired_counts["VH169"])
current_counts <- table(factor(combined_intensities$classification, levels = names(desired_counts)))
if (current_counts["564"] < min_564) {
  deficit <- min_564 - current_counts["564"]
  
  # Select closest cells from WT and VH169 to move to 564
  candidates <- combined_intensities %>%
    filter(classification %in% c("WT", "VH169")) %>%
    mutate(
      distance_to_564 = abs(cd45_1_intensity - thresholds$`564`$cd45_1_intensity) +
                        abs(cd45_2_intensity - thresholds$`564`$cd45_2_intensity)
    ) %>%
    arrange(distance_to_564)
  
  # Reassign candidates iteratively
  for (i in seq_len(nrow(candidates))) {
    candidate <- candidates[i, ]
    candidate_class <- candidate$classification
    
    # Check if moving this candidate would violate WT or VH169
    if (candidate_class == "WT" && current_counts["WT"] - 1 < min_WT) next
    if (candidate_class == "VH169" && current_counts["VH169"] - 1 < min_VH169) next
    
    # Reassign candidate to 564
    combined_intensities$classification[combined_intensities$label == candidate$label] <- "564"
    current_counts["564"] <- current_counts["564"] + 1
    current_counts[candidate_class] <- current_counts[candidate_class] - 1
    
    # Stop once deficit is resolved
    if (current_counts["564"] >= min_564) break
  }
}

# Print final percentages
final_percentages <- combined_intensities %>%
  count(classification) %>%
  mutate(percentage = n / sum(n) * 100)
print(final_percentages)

# Define annotation colors
annotation_mapping_colors <- list(
  "WT" = c(144, 238, 144),        # Light green
  "564" = c(0, 0, 139),           # Dark blue
  "F_VH169" = c(139, 0, 0),       # Dark red
  "L_VH169" = c(255, 182, 193),   # Light red
  "Unclassified" = c(0, 0, 0)     # Black
)

# Annotate each sample and generate annotated CSV, cell counts, and corresponding images
for (prefix in sample_prefixes) {
  CD45_1_file <- paste0(prefix, "_CD45_1.csv")
  CD45_2_file <- paste0(prefix, "_CD45_2.csv")
  G6_file <- paste0(prefix, "_G6.csv")
  mask_file <- paste0(prefix, "_cp_masks.png")
  if (all(c(CD45_1_file, CD45_2_file, G6_file, mask_file) %in% files)) {
    output_annotated_csv <- paste0(prefix, "_annotated_cells.csv")
    output_counts <- paste0(prefix, "_annotated_counts.csv")
    output_image <- paste0(prefix, "_annotated_image.png")
    
    # Load region properties
    CD45_1_props <- read_csv(CD45_1_file, show_col_types = FALSE)
    CD45_2_props <- read_csv(CD45_2_file, show_col_types = FALSE)
    G6_props <- read_csv(G6_file, show_col_types = FALSE)
    
    # Ensure unique labels by prefixing them
    CD45_1_props$label <- paste0(prefix, "_", CD45_1_props$label)
    CD45_2_props$label <- paste0(prefix, "_", CD45_2_props$label)
    G6_props$label <- paste0(prefix, "_", G6_props$label)
    
    # Rename G6 mean_intensity to G6_intensity for consistency
    G6_props <- G6_props %>%
      rename(G6_intensity = mean_intensity)
    
    # Merge the tables by label
    merged_props <- CD45_1_props %>%
      inner_join(CD45_2_props, by = "label", suffix = c("_CD45_1", "_CD45_2")) %>%
      inner_join(G6_props, by = "label")
    
    # Assign annotations from global classification
    merged_props <- merged_props %>%
      mutate(annotation = combined_intensities$classification[match(label, combined_intensities$label)])
    
    # Distinguish F vs L based on G6 intensity
    merged_props <- merged_props %>%
      mutate(
        F_L = ifelse(
          annotation %in% c("564", "WT"),  # Check if the cell type is 564 or B6
          annotation,  # Keep the original annotation for 564 and B6
          paste0(
            ifelse(G6_intensity >= threshold_G6, "F_", "L_"), 
            annotation  # Add prefix only for non-564 and non-B6 cells
          )
        )
      )
    
    # Save the annotated data as CSV
    write_csv(merged_props, output_annotated_csv)
    
    # Count the number of cells for each annotation
    cell_type_counts <- merged_props %>%
      group_by(F_L) %>%
      summarise(cell_count = n()) %>%
      ungroup()
    
    # Save the counts as CSV
    write_csv(cell_type_counts, output_counts)
    
    # Load the corresponding segmentation mask
    segmentation_mask <- load.image(mask_file)
    segmentation_matrix <- as.matrix(segmentation_mask)
    
    # Scale up values in the segmentation matrix for 16-bit intensity range
    segmentation_matrix <- segmentation_matrix * 65535
    
    # Initialize RGB channels for the output image
    red_channel <- matrix(0, nrow = nrow(segmentation_matrix), ncol = ncol(segmentation_matrix))
    green_channel <- matrix(0, nrow = nrow(segmentation_matrix), ncol = ncol(segmentation_matrix))
    blue_channel <- matrix(0, nrow = nrow(segmentation_matrix), ncol = ncol(segmentation_matrix))
    
    # Assign colors to the matrix based on annotation
    for (i in 1:nrow(merged_props)) {
      label <- merged_props$label[i]
      annotation <- merged_props$F_L[i]
      if (!is.null(annotation_mapping_colors[[annotation]])) {
        
        # Extract the numeric part of the label for matching
        numeric_label <- as.numeric(sub(".*_", "", label))
        
        # Find matching indices in the segmentation matrix
        matching_indices <- which(segmentation_matrix == numeric_label)
        if (length(matching_indices) > 0) {
          color <- annotation_mapping_colors[[annotation]]
          red_channel[matching_indices] <- color[1]
          green_channel[matching_indices] <- color[2]
          blue_channel[matching_indices] <- color[3]
        }
      } else {
        print(paste("No mapping found for annotation:", annotation))
      }
    }
    
    # Combine RGB channels into a single 3D image and scale the values between 0 and 1
    annotated_image <- array(0, dim = c(nrow(segmentation_matrix), ncol(segmentation_matrix), 3))
    annotated_image[,,1] <- red_channel / 255
    annotated_image[,,2] <- green_channel / 255
    annotated_image[,,3] <- blue_channel / 255
    
    # Flip the image horizontally
    annotated_image <- annotated_image[, ncol(annotated_image):1, ]
    
    # Rotate the annotated image 90 degrees counterclockwise
    annotated_image <- aperm(annotated_image, c(2, 1, 3))
    annotated_image <- annotated_image[nrow(annotated_image):1, , ]
    
    # Save the annotated image as PNG
    writePNG(annotated_image, output_image)
    cat("Annotated image saved for", prefix, "\n")
    
    # Clear variables to free memory
    rm(CD45_1_props, CD45_2_props, G6_props, merged_props, segmentation_mask,
       segmentation_matrix, annotated_image, cell_type_counts, red_channel,
       green_channel, blue_channel)
  }
}

# Run garbage collection
invisible(gc())
print("Compeleted!")

