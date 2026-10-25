Upper Vaal Water Quality Spatial Statistics Dashboard
Project Overview
This project presents an interactive R Shiny dashboard for analysing water quality in the Upper Vaal Water Management Area, South Africa, from 2016 to 2020.
The dashboard uses water quality observations from 45 monitoring stations to investigate spatial patterns and changes over time.
Water Quality Parameters
pH
Electrical Conductivity (EC)
Sulphate (SO4)
Annual median values are used to represent each parameter at each monitoring station.
Spatial Statistics Methods
1. Global Moran's I
Measures whether water quality values show spatial clustering, dispersion or a random spatial pattern across the study area.
2. Local Moran's I (LISA)
Identifies local clusters and spatial outliers, including High–High, Low–Low, High–Low and Low–High patterns.
3. Ordinary Kriging
Uses spatial interpolation to estimate water quality values at locations without monitoring stations.
Dashboard Features
Interactive maps of the Upper Vaal Water Management Area
Year selection from 2016 to 2020
Water quality parameter selection
Global Moran's I results
Local Moran's I cluster maps
Ordinary Kriging interpolation maps
Monitoring station information
Water quality trends
Additional GIS layers, including catchments, rivers, dams and land use
Technology
Developed using R, RStudio, Shiny, sf, leaflet, dplyr, ggplot2 and other supporting R packages.
Live Dashboard
https://spatialstats.shinyapps.io/upper-vaal-water-quality/
Academic Purpose
Developed as part of a Spatial Statistics group assignment at the University of the Witwatersrand.
 