# Herbivore population differences rival geographic and biophysical variation in structuring ecosystem function
*Global Change Biology* (2025)

Authors: Nathalie R. Sommer, Annise M. Dobson, Matthew S. Baker, Geoffrey C. Trussell, and Oswald J. Schmitz

Code by Nathalie Sommer

This repository contains all data, scripts, and documentation necessary to reproduce the analyses and figures presented in our manuscript. Upon acceptance, this README will be updated with the manuscript DOI.

## Table of Contents
- [Project Overview](#project-overview)
- [Data](#data)
- [Scripts](#scripts)
- [Supporting Information](#supporting-information)
- [Dependencies](#dependencies)
- [Reproducibility](#reproducibility)

## Project Overview
Geographic variation in ecosystem function is often attributed to differences in climate and soil properties, with biophysical constraints assumed to dictate spatial patterns in nutrient cycling, carbon storage, and plant productivity. However, biotic interactions, particularly herbivory, also vary geographically and can generate feedbacks that influence ecosystem processes. Using a replicated three-year field experiment, we tested how population-level functional differences in a widespread arthropod herbivore mediate geographic variation in ecosystem function. Structural equation modeling revealed that herbivores exerted strong direct effects on plant biomass, soil carbon, and nitrogen mineralization, often surpassing the influence of historical conditions and geographic variation in climate. Moreover, functionally distinct herbivore populations had divergent effects on nutrient cycling and plant diversity, demonstrating that population-level differences introduce novel pathways of influence on ecosystem function. These findings challenge ecosystem models that prioritize abiotic constraints and highlight the need to incorporate consumer-driven feedbacks into ecological frameworks.


## Data
All data is contained within the `Data` folder in their raw, original form. The data is organized into the following subdirectories:

### Soil Properties
- `BulkDensity/`: Soil bulk density measurements
- `GravimetricWater/`: Soil moisture content data
- `LOI/`: Loss on ignition data for organic matter content
- `N-min/`: Nitrogen mineralization data
- `SIR/`: Substrate-induced respiration measurements for microbial biomass
- `Texture_pH/`: Soil texture and pH measurements
- `WaterHoldingCapacity/`: Soil water holding capacity data

### Vegetation Data
- `VegBiomass/`: Plant biomass measurements
- `VegCommunity/`: Plant community composition data

### Chemical Analysis
- `Elements/`: Elemental analysis data

## Supporting Information
This folder contains the script for the accompanying [Shiny App](https://nathaliesommer.shinyapps.io/herbivore-populations-structure-ecosystems/) and a script for a multigroup SEM analysis. 

## Scripts
All analysis scripts are contained within the `Scripts` folder:

### Data Processing and Analysis
- `LOI.R`: Processes loss on ignition data
- `N-min.R`: Analyzes nitrogen mineralization rates
- `SIR.R`: Processes substrate-induced respiration data
- `SoilPlantCN.R`: Analyzes soil and plant C:N data
- `Texture-BD-pH.R`: Processes soil physical properties
- `Veg.R`: Analyzes vegetation data
- `SEM_Composite.R`: Main structural equation modeling analysis. The initial SEM results are stored in `SEM_Summary_Initial.txt`, which contains the first iteration of the structural equation model prior to integrating the tests of directed separation.

## Dependencies
The analyses were conducted in R (version 4.4.1) with the following packages:

### Data Cleaning and Analysis
- `tidyverse` (v1.3.0): Collection of data manipulation packages
- `dplyr` (v1.1.4): Data manipulation
- `tidyr` (v1.3.0): Data tidying
- `broom` (v1.0.5): Converting statistical objects into tidy data frames
- `purrr` (v1.0.2): Functional programming tools

### Statistical Analysis
- `lme4` (v1.1-35.5): Linear mixed-effects models
- `DHARMa` (v0.4.6): Residual diagnostics for hierarchical models
- `car` (v3.1-30): Companion to Applied Regression
- `boot` (v1.3-30): Bootstrap methods
- `piecewiseSEM` (v2.3.0): Piecewise structural equation modeling
- `MuMIn` (v1.47.5): Multi-model inference

### Visualization
- `ggplot2` (v3.5.1): Data visualization
- `ggforce` (v0.4.1): Extensions to ggplot2
- `ggdag` (v0.2.10): Visualization of directed acyclic graphs
- `DiagrammeR` (v1.0.11): Graph and network visualization
- `RColorBrewer` (v1.1-3): Color palettes
- `paletteer` (v1.5.0): Comprehensive collection of color palettes

## Reproducibility
To reproduce the analysis and figures:

1. Clone the repository:
   ```bash
   git clone [repository URL]
   ```

2. Set your working directory to the cloned repository

3. Install required R packages:
   ```R
   # Code for package installation will be provided
   ```

4. Run the scripts in the following order:
   - First run data processing scripts:
     - `LOI.R`
     - `N-min.R`
     - `SIR.R`
     - `Texture-BD-pH.R`
   - Then run the main analysis:
     - `SEM_Composite.R`

Results will be saved in the appropriate output directories as specified in each script.

## Contact
For questions about the code or data, please open an issue in this repository or contact Nathalie Sommer.
