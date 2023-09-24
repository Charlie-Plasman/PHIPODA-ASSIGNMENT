# PHIPODA-ASSIGNMENT
Quality control, visualisations and insights on biodiversity data for a job interview assignment.

# README

## Introduction

This README provides an overview of the R script and data processing steps used for the assignment related to the "SCIENTIFIC COLLABORATOR ANTARCTIC AND SOUTHERN OCEAN BIODIVERSITY OPPORTUNITY" offered by the RBINS. The script performs various data preprocessing, validation, and cleaning.

## Prerequisites

Before running the R script, ensure you have the following installed:

- R (https://www.r-project.org/)
- Required R packages (`readr`, `dplyr`, `SOmap`, `raadtools`, `raster`, and `worrms`). You can install these packages using the `install.packages()` function.

## Usage

1. **Data Preparation**: 

   - Place your dataset file (`occurrence.csv`, available in this repository) in a chosen directory.
   - Update the file path in the script.

2. **Data Preprocessing**:

   - The script loads the dataset, checks data types, and identifies missing data.
   - It validates latitude, longitude, and depth values and filters out invalid records.
   - Potential duplicate records are removed based on the "ID" column.
   - The script categorizes sampling methods and checks for consistency in the dynamic properties variable.

3. **Spatial Filtering**:

   - The code applies spatial filters to remove data points on land using a land mask.
   - It also checks for inconsistencies in dynamic properties and filters data based on keywords.
    
4. **WoRMS Matching**:

   - The code uses the 'worrms' package to match scientific names with WoRMS.
   - It assigns WoRMS identifiers (aphiaID) and valid names to the dataset.
   - The script checks for multiple matches and prioritizes accepted names.

5. **Results**:

   - The script provides summary statistics, frequency counts of sampling method and several biodiversity data visualizations.

## Output

- The processed dataset is cleaned for missing/inconsistent data, updated with WoRMS identifiers/valid names and aggregated sampling methods.
- Maps for amphipoda biodiversity data visualizations are also available in the repository.
- Results are displayed in the R console, and you can customize further actions as needed.
  

## Notes

- Some data preprocessing steps may require manual review or further WoRMS corrections.
- Make sure to adapt the script according to your specific dataset and requirements.

## License

This script is provided under the [MIT License](LICENSE).

## Author

Charlie Plasman

## Contact

For any questions or issues, please contact charlie.plasman@umb.be.

---
