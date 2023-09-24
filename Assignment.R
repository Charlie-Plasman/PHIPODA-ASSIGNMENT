# PLASMAN CHARLIE
# 09/2023
################################################################################



### Accès aux données
library(readr)
setwd("~/Desktop/RBINS Job/Assignment")
data <- read_delim("occurrence.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(data)
# Check data types
str(data)
# Check for missing data
missing_data <- colSums(is.na(data))
missing_data <- missing_data[missing_data > 0]
print(missing_data)



### Filtres spatiaux
# On arrondit à 4 décimales, pourquoi ?
data$decimalLatitude = round(data$decimalLatitude, 4)
data$decimalLongitude = round(data$decimalLongitude, 4)
# Coordinate validation (assuming latitude and longitude columns are named "decimalLatitude" and "decimalLongitude")
# Check for out-of-range latitude values
invalid_lat <- subset(data, decimalLatitude < -90 | decimalLatitude > -35)
print(invalid_lat)
# Check for out-of-range longitude values
invalid_lon <- subset(data, decimalLongitude < -180 | decimalLongitude > 180)
print(invalid_lon)
# Depth data validation (assuming "Minimum Depth in Meters" and "Maximum Depth in Meters" columns)
# Check if minimum depth is greater than maximum depth
invalid_depth <- subset(data, `minimumDepthInMeters` > `maximumDepthInMeters`)
print(invalid_depth)
library(dplyr)
data <- data %>%
  filter(!(data$id %in% invalid_depth$id))
# Filter out invalid observations from data
data <- data %>%
  filter(!(data$id %in% invalid_lat$id))
# Vérifier pour le point au lat et long égales
equivalent = (data$decimalLatitude == data$decimalLongitude)
equi_points = data[equivalent,]


library(SOmap)
basemap = SOmap(bathy_legend = F, ice = T, trim = -35, fronts = "Park", border = F) 
# Load land mask data
land_mask <- raster("Antarctica_land.nc")
# Extract longitude and latitude values
spat_occ = cbind.data.frame(data$decimalLongitude,data$decimalLatitude)
colnames(spat_occ) = c("longitude", "latitude")
head(spat_occ)
plot(basemap)
SOplot(spat_occ, col = "darkblue", cex = 0.5, pch =20)
# Create a logical vector indicating if each point is at sea
library(raadtools)
on_land <- !is.na(extract(land_mask, test))
data <- filter(data, !on_land)
plot(basemap)
SOplot(cbind.data.frame(data$decimalLongitude,data$decimalLatitude), col = "darkblue", cex = 0.5, pch =20)
# On passe de 1808 points à 1776 avec le masque terrestre




### Filtres temporels
# Regarder si moyen de retrouver l'info temporelle pour le sampling de chaque ligne

### Filtres taxonomiques
{
library(tidyverse)
library(worrms)

taxon_match_fn <- function(unames) {
  
  pt <- proc.time()
  
  # split by blocks of 50
  pages <- split(unames, as.integer((seq_along(unames) - 1) / 50))
  
  # for each block of 50, use worrms
  taxon_match <- pages %>% 
    map_df(., function(unames_i) {
      
      tryCatch({ 
        result <- worrms::wm_records_taxamatch(unames_i) # Get records for one or more taxonomic name(s) using the TAXAMATCH fuzzy matching algorithm
      }, error=function(e){result<-NULL}) #print(paste(p, unames_i)); 
      
      if (!is.null(result)) {
        result_df <- map_df(1:length(result), function(i) {
          if (nrow (result[[i]]) > 0) {
            df <- result[[i]]
            df$input_name <- unames_i[i] # associate the input name to each match (because sometimes what is returned in the scientificname column is not exactly the input name)
            df$nb_match <- nrow(df) # to identify multiple matches
            df$nb_accepted_match <- sum(df$status == "accepted") #to identify multiple matches left after removing unaccepted/uncertain names
            return(df)
          }
        })
      } else {result_df = NULL}
      return(result_df)
    })
  
  print(proc.time() - pt)
  gc() 
  
  return(taxon_match)
}

taxon_match_iter <- function(all_names) {
  n_taxa <- length(all_names)
  pt_match <- proc.time()
  for (m in 1:ceiling(n_taxa/block_size)) {
    start <- (m-1)*block_size+1
    end <- min(c(m*block_size, n_taxa))
    tax_mtch_sub <- taxon_match_fn(all_names[start:end])
    if (nrow(tax_mtch_all) == 0) {
      tax_mtch_all <- tax_mtch_sub
    }
    tax_mtch_all <- rbind(tax_mtch_all, tax_mtch_sub)
    print(m*block_size)
    print(proc.time() - pt_match)
  }
  return(tax_mtch_all)
}

tax_match_from_file <- function(input_file, output_file, block_size = 2000) {
  all_names <- read.csv2(input_file)$scientificName %>% 
    as.character()
  
  tax_mtch_all <- NULL
  tax_mtch_all <- taxon_match_iter(all_names)
  
  tax_mtch_all %>% write_csv(output_file)
  
  for (k in 1:10) {
    failed_attempts <- setdiff(all_names, tax_mtch_all$input_name)
    length(failed_attempts)
    if (length(failed_attempts) > 0) { # identify the ones to retry in case we lost the internet connection along the way
      tax_mtch_all <- taxon_match_iter(failed_attempts)
    } 
  }
  tax_mtch_all %>% write_csv(output_file) 
}

block_size = 2000
}
# Call the tax_match_from_file function
tax_match_from_file("occurrence.csv", "output_matches.csv", block_size)


# Initialize an empty vector to store aphiaIDs
aphiaIDs <- vector("numeric", length = nrow(data))

# Loop through each row of your dataset and match with WoRMS
for (i in 1:nrow(data)) {
  scientific_name <- data$scientificName[i]
  
  # Use the worrms_lookup() function to match the name with WoRMS
  result <- wm_records_taxamatch(name = scientific_name)
  
  # Check if a match was found
  if (!is.null(result$aphiaID)) {
    aphiaIDs[i] <- result$aphiaID
  } else {
    # Handle unmatched names here (you can flag them for further review)
    aphiaIDs[i] <- NA
  }
}

# Add the aphiaIDs to your dataset
data$aphiaID <- aphiaIDs

# Print the first few rows of your updated dataset
head(data)


### Filtres "types d'échantillonnage"
# Check for inconsistencies in dynamic properties
# Check for spelling and formatting consistency
unique_methods <- data.frame(unique(data$`dynamicProperties`))
print(unique_methods)
# Check for missing or unspecified methods
missing_methods <- subset(data, is.na(dynamicProperties) | dynamicProperties == "" | dynamicProperties == "	
{gear:}")
print(missing_methods)
# Filter out missing methods from data
data <- data %>%
  filter(!(data$id %in% missing_methods$id))

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
data <- data %>%
  mutate(methodCategory = ifelse(grepl("trawl", dynamicProperties, ignore.case = TRUE), "Trawl",
                         ifelse(grepl("trap", dynamicProperties, ignore.case = TRUE), "Trap",
                                ifelse(grepl("scuba", dynamicProperties, ignore.case = TRUE), "Scuba",
                                       ifelse(grepl("net", dynamicProperties, ignore.case = TRUE), "Net",
                                              ifelse(grepl("dredge", dynamicProperties, ignore.case = TRUE), "Dredge", "Others"))))))
# Avec la fonction, j'ai rassemeblé en 10 catégories et ensuite réduit manuellement à 5 (plus de 30 occurrences)
# Frequency of methods
method_counts <- table(data$methodCategory)
print(method_counts)

### Filtres de duplication
# Remove duplicate records based on the ID column (assuming "ID" is the unique identifier)
data <- data[!duplicated(data$id), ]
identical(data$id, data$occurrenceID)
duplicate_rows <- duplicated(data[c("decimalLatitude", "decimalLongitude","scientificName")])
duplicated_data = data[duplicate_rows,]
# ATTENTION, certains duplicats sont seuls ce qui n'a pas de sens et regarder après correction worms parce que certains noms
