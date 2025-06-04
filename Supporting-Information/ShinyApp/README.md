# Herbivore Populations Restructure Ecosystems

Interactive visualization app for exploring the effects of herbivore populations on ecosystem structure, based on research from Sommer et al. (Global Change Biology).

## Overview

This Shiny application provides interactive visualizations of how different herbivore populations (Reactor, Resistor, and Vegetation treatments) affect various ecosystem metrics across different sites and years.

## Features

### Raw Data Visualization
- Explore relationships between any two ecosystem metrics
- Filter data by treatment type, year, and site
- View regression lines with confidence intervals
- Compare patterns across different years

### Common Garden Comparison
- Compare treatment effects across different sites and years
- View statistical summaries and ANOVA results
- Analyze treatment effects on various ecosystem metrics

## Available Metrics

The app includes the following ecosystem metrics:
- Goldenrod Biomass (g)
- Grass Biomass (g)
- Plant Diversity (Shannon-Wiener)
- Soil %C and %N
- SIR (CO2 per hour per g)
- N-mineralization (mg N/cm² per month)
- Litter %C and %N
- Plant tissue %C and %N

## Running the App Locally

1. Ensure you have R installed (version 4.0.0 or higher recommended)
2. Install required packages:
   ```R
   install.packages(c("shiny", "dplyr", "ggplot2", "DT", "readr", "here", 
                     "lme4", "lmerTest", "emmeans", "multcomp"))
   ```
3. Clone this repository
4. Open `app.R` in RStudio
5. Click "Run App" or use `shiny::runApp()`

## Data Sources

The app uses processed data from the following sources:
- Vegetation community composition
- Soil and plant C/N content
- Soil respiration (SIR)
- N-mineralization rates
- Soil physical properties

## Citation

If you use this app in your research, please cite:
Sommer et al. Herbivore population differences rival geographic and biophysical variation in structuring ecosystem function. Global Change Biology.

## Contact

For questions or issues, please contact Nathalie Sommer (nathalie.sommer@aya.yale.edu)