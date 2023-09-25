################################################################################
# PLASMAN CHARLIE - SEPTEMBER 2023
# SCIENTIFIC COLLABORATOR ANTARCTIC AND SOUTHERN OCEAN BIODIVERSITY 
# RBINS JOB OPPORTUNITY AND ASSIGNMENT
################################################################################



################################################################################ Data access and overview
# Load libraries
library(readr)
library(dplyr)
library(SOmap)
library(raadtools)
library(raster)
library(worrms)
library(ggplot2)

# Set the working directory
setwd("~/Desktop/RBINS Job/Assignment")

# Read data from CSV file
data <- read_delim("occurrence.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

# View the data
View(data)

# Check for missing data
missing_data <- colSums(is.na(data))
missing_data <- missing_data[missing_data > 0]
print(missing_data)

################################################################################ Spatial Filters
# Round latitude and longitude to 4 decimal places
data$decimalLatitude <- round(data$decimalLatitude, 4)
data$decimalLongitude <- round(data$decimalLongitude, 4)

# Validate latitude and longitude values
invalid_lat <- subset(data, decimalLatitude < -90 | decimalLatitude > -35)
invalid_lon <- subset(data, decimalLongitude < -180 | decimalLongitude > 180)
print(invalid_lat)
print(invalid_lon)

# Depth data validation
invalid_depth <- subset(data, minimumDepthInMeters > maximumDepthInMeters)
print(invalid_depth)

# Filter out invalid observations based on depth
data <- data %>%
  filter(!(data$id %in% invalid_depth$id))

# Filter out invalid observations based on latitude
data <- data %>%
  filter(!(data$id %in% invalid_lat$id))

# Filter out equivalent points where latitude equals longitude
equivalent <- (data$decimalLatitude == data$decimalLongitude)
equi_points <- data[equivalent,]
data <- data %>%
  filter(!(data$id %in% equi_points$id))

# Create a basemap
basemap <- SOmap(bathy_legend = FALSE, ice = TRUE, trim = -35, fronts = "Park", border = FALSE) 

# Load land mask data
land_mask <- raster("Antarctica_land.nc")

# Extract longitude and latitude values for spatial data
spat_occ <- cbind.data.frame(data$decimalLongitude, data$decimalLatitude)
colnames(spat_occ) <- c("longitude", "latitude")

# Plot the basemap
plot(basemap)

# Plot spatial points on the basemap
SOplot(spat_occ, col = "darkblue", cex = 1, pch = 20)

# Create a logical vector indicating if each point is at sea
on_land <- !is.na(extract(land_mask, spat_occ, method = 'simple'))

# Filter out points on land
data <- filter(data, !on_land)

# Plot the updated basemap
plot(basemap)

# Plot updated spatial points on the basemap
SOplot(cbind.data.frame(data$decimalLongitude, data$decimalLatitude), col = "darkblue", cex = 1, pch = 20)


################################################################################ Taxonomic Filters
# Alternative 1 : Fill records missing a scientificNameID and check for preferred accepted taxon name for all records
# Define a function to match scientific names with WoRMS for a subset of the data
match_with_worms_subset <- function(dataset, start_row, end_row) {
  results_list <- list()
  aphiaIDs <- vector("numeric", length = nrow(dataset))
  valid_names <- character(length = nrow(dataset))
  
  for (row_index in start_row:end_row) {
    scientific_name <- dataset$scientificName[row_index]
    result <- wm_records_taxamatch(scientific_name)
    results_list[[row_index]] <- result
    
    for (record_index in seq_along(result)) {
      if (!is.null(result[[record_index]]$AphiaID)) {
        accepted_index <- which(result[[record_index]]$status == "accepted")
        if (length(accepted_index) > 0) {
          aphiaIDs[row_index] <- result[[record_index]]$AphiaID[accepted_index[1]]
          valid_names[row_index] <- result[[record_index]]$valid_name[accepted_index[1]]
        } else {
          aphiaIDs[row_index] <- result[[record_index]]$AphiaID[1]
          valid_names[row_index] <- result[[record_index]]$valid_name[1]
        }
      } else {
        aphiaIDs[row_index] <- NA
        valid_names[row_index] <- NA
      }
    }
  }
  return(list(results_list = results_list, aphiaIDs = aphiaIDs, valid_names = valid_names))
}

# Get records with missing scientificNameID
missing_IDs <- subset(data, is.na(data$scientificNameID))

# Define the size of the data processing block (e.g., 100 rows per block)
block_size <- 67

# Initialize empty vectors for aphiaIDs and valid names
aphiaIDs <- vector("numeric", length = nrow(missing_IDs))
valid_names <- character(length = nrow(missing_IDs))

# Process the data in blocks
for (start_row in seq(1, nrow(data), by = block_size)) {
  end_row <- min(start_row + block_size - 1, nrow(data))
  
  # Match scientific names with WoRMS for the current block
  matching_results <- match_with_worms_subset(missing_IDs, start_row, end_row)
  
  # Update the aphiaIDs and valid names vectors with results from the current block
  aphiaIDs[start_row:end_row] <- matching_results$aphiaIDs[start_row:end_row]
  valid_names[start_row:end_row] <- matching_results$valid_names[start_row:end_row]
}

# Add the aphiaIDs and valid names to your dataset
missing_IDs$scientificNameID <- aphiaIDs
missing_IDs$scientificName <- valid_names

# Print the first few rows of your updated dataset
head(missing_IDs)

# Alternative 2 : Remove records with missing scientificNameID
# Filter out observations without a taxon name ID
# Get records with missing scientificNameID
missing_IDs <- subset(data, is.na(data$scientificNameID))
data = data <- data %>%
  filter(!(data$id %in% missing_IDs$id))

################################################################################ Observation type filters
# Check for inconsistencies in dynamic properties (sampling methods)
unique_methods <- unique(data$dynamicProperties)
print(unique_methods)

# Check for missing or unspecified methods
missing_methods <- subset(data, is.na(dynamicProperties) | dynamicProperties == "" | dynamicProperties == "{gear:}")
print(missing_methods)

# Filter out missing methods from data
data <- data %>%
  filter(!(id %in% missing_methods$id))

# Divide the data into "scuba diving" and "non-scuba" categories
data$observation_type <- ifelse(grepl("scuba", data$dynamicProperties, ignore.case = TRUE), "scuba diving", "bottom trawling")

# Frequency of methods
method_counts <- table(data$observation_type)
print(method_counts)

# Bonus : Quick analysis of the sampling methods
# Define a function to filter data based on a keyword
filter_by_keyword <- function(data, keyword) {
  data %>%
    filter(grepl(keyword, dynamicProperties, ignore.case = TRUE))
}

# List of keywords
keywords <- c("trawl", "grab", "trap|traps", "scuba", "dredge", "corer", "kelp", "sledge", "fish", "net")

# Create a list of filtered data frames
filtered_data <- lapply(keywords, function(keyword) {
  filter_by_keyword(data, keyword)
})

# Categorize methods
data <- data %>%
  mutate(
    methodCategory = ifelse(grepl("trawl", dynamicProperties, ignore.case = TRUE), "Trawl",
                            ifelse(grepl("trap", dynamicProperties, ignore.case = TRUE), "Trap",
                                   ifelse(grepl("scuba", dynamicProperties, ignore.case = TRUE), "Scuba",
                                          ifelse(grepl("net", dynamicProperties, ignore.case = TRUE), "Net",
                                                 ifelse(grepl("dredge", dynamicProperties, ignore.case = TRUE), "Dredge", "Others"))))))

# Frequency of methods
method_counts <- table(data$methodCategory)
print(method_counts)

################################################################################ Duplication filters
# Check for duplicates in "id" and equivalence between "id" and "occurenceID"
data <- data[!duplicated(data$id), ]
identical(data$id, data$occurrenceID)

# Identify duplicates on all rows 
duplicate_rows <- duplicated(sorted_data[c("decimalLatitude", "decimalLongitude")])

# Filter data to keep only duplicates on all selected columns 
all_columns <- c("decimalLatitude", "decimalLongitude", "scientificName","scientificNameID","minimumDepthInMeters","maximumDepthInMeters")
all_duplicate_rows <- duplicated(sorted_data[all_columns])
duplicated_data <- sorted_data[all_duplicate_rows, ]

# Filter out potential duplicated data
data <- data %>%
  filter(!(data$id %in% duplicated_data$id))

################################################################################  Processed data save
# Save the processed data as a .csv file
write.csv(data, "processed_data.csv", row.names = FALSE)

################################################################################  Visualization
# Plot the updated basemap
plot(basemap)

# Plot updated spatial points on the basemap
SOplot(cbind.data.frame(data$decimalLongitude, data$decimalLatitude), col = "darkblue", cex = 1, pch = 20)

# Create a logical vector indicating the observation types 
is_scuba_diving <- data$observation_type == "scuba diving"

# Plot a map with the two types of observations highlighted
plot(basemap)
SOplot(data[is_scuba_diving, c("decimalLongitude", "decimalLatitude")], col = "darkorange", cex = 1, pch = 20)
SOplot(data[!is_scuba_diving, c("decimalLongitude", "decimalLatitude")], col = "darkgreen", cex = 1, pch = 20)

# Create a zoomed map with the two types of observations highlighted
depth = raster("depth.nc")
extent = extent(-90, 0, -90, -25)
zoom = crop(depth, extent)

# Divide pair of coordinates in two groups
trawling = data[!is_scuba_diving, c("decimalLongitude", "decimalLatitude")]
diving = data[is_scuba_diving, c("decimalLongitude", "decimalLatitude")]

# Plot observations with two colors
plot(zoom, col = hcl.colors(80, "Blues"), legend = F)
plot(SpatialPoints(trawling), col = "red", cex = 2, pch = 20, add = T)
plot(SpatialPoints(diving), col = "blue", cex = 2, pch = 20, add = T)

# Create a data frame with depth ranges
depth_data <- data.frame(
  observation = 1:nrow(data),  # Unique identifier for each observation
  min_depth = -data$minimumDepthInMeters,
  max_depth = -data$maximumDepthInMeters,
  observation_type = data$observation_type 
)

# Calculate mean depth
mean_depth <- mean(depth_data$min_depth + depth_data$max_depth, na.rm = T) / 2
upper_depth <- quantile(depth_data$min_depth + depth_data$max_depth, probs = 0.90, na.rm = TRUE) / 2
lower_depth <- quantile(depth_data$min_depth + depth_data$max_depth, probs = 0.10, na.rm = TRUE) / 2

# Create a ggplot plot for depth ranges 
ggplot(depth_data, aes(x = observation, ymin = min_depth, ymax = max_depth, color = observation_type)) +
  geom_linerange(size = 1) +
  geom_hline(yintercept = mean_depth, linetype = "dashed", color = "darkgrey", size = 1) +  # Add mean line
  geom_hline(yintercept = upper_depth, linetype = "dashed", color = "black", size = 1) +  # Add upper line
  geom_hline(yintercept = lower_depth, linetype = "dashed", color = "black", size = 1) +  # Add lower line
  labs(title = "Depth range (processed dataset)",
       x = "Observations",
       y = "Depth (meters)") +
  scale_color_manual(values = c("scuba diving" = "blue", "bottom trawling" = "red"), name = "Observation types") +
  theme_minimal() +
  theme(axis.text.x = element_blank())  # Hide x-axis labels

# Create a KDE layer (proxy sampling effort)
KDE_layer = raster(MASS::kde2d(data$decimalLongitude, data$decimalLatitude,
                               n=c(ncol(depth), nrow(depth)),
                               lims=c( -180, 180, -90, -35))) 
extent(KDE_layer) = extent(depth)
KDE_layer = mask(KDE_layer, depth) 

# Plot the KDE mask onto the zoomed map
plot(zoom, legend = F, col = "lightblue")
plot(KDE_layer, legend = T, add = T, legend) 
