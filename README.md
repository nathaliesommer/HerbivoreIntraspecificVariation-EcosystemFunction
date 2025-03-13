# Herbivore population differences rival geographic and biophysical variation in structuring ecosystem function
[Manuscript in review]

Authors: Nathalie R. Sommer, Annise M. Dobson, Matthew S. Baker, Geoffrey C. Trussell, and Oswald J. Schmitz

Code by Nathalie Sommer

This repository contains all data, scripts, and documentation necessary to reproduce the analyses and figures presented in our manuscript. Upon acceptance, this README will be updated with the manuscript DOI.

## Table of Contents
- [Project Overview](#project-overview)
- [Data](#data)
- [Scripts](#scripts)
- [Dependencies](#dependencies)
- [Reproducibility](#reproducibility)

## Project Overview
We report on findings from a three-year geographically replicated experiment designed to quantify how differences in trait expressions of local populations of a dominant herbivore species mediate variation in ecosystem function across its geographic range in the New England region of the eastern USA. The species—Melanoplus femurrubrum—is a widespread grasshopper in which local populations differently express physiological and behavioral traits (Parsons & Joern 2014; Rosenblatt et al. 2019; Baker et al; Sommer et al.). We applied structural equation modeling (SEM) to quantify the relative contributions of among-population variation, and local climate and historical legacies as drivers of key ecosystem variables and functions, including plant biomass, soil nutrients, and nitrogen mineralization. We show how geographic variation in grasshopper impacts on old field ecosystem functioning varies with population origin (warm vs cool sites) which differ systematically in daily mean maximum temperature during the growing season (Baker et al.). Differences in old field ecosystem impacts between replicated warm and cool site populations can be attributed to differences in their expression of grasshopper population physiological and behavioral trait plasticity under environmental stressors. This study shows how mechanisms determining within-species variation in herbivory can manifest spatially to play an important role in shaping geographic variation in ecosystem function.

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

### Data Manipulation and Analysis
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
For questions about the code or data, please open an issue in this repository or contact [contact information].

## License
[License information to be added]