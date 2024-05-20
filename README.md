# PCR Analyzer & Visualizer

## Overview
PCR Analyzer & Visualizer is a  R application designed to facilitate the analysis and visualization of Polymerase Chain Reaction (PCR) normalized results.
## Features
- **Data Import**: Easily import PCR data from Excel or CSV. 
- **Analysis**: Conduct analysis of PCR results. The output is similar to the easyqPCR's one - a list with expression values and sd.
- **Visualization**: Using ggplot2 and plotly, generate high-quality plots and graphs to visually represent PCR data.
- **Export**: Export visualizations in multiple formats for further analysis or presentation.

## Installation
1. Clone or download this repository to your local machine.
2. Ensure you have R and RStudio installed on your system.
3. Open `PCR_Analyzer_Visualizer.Rproj` in RStudio.
4. Install required packages by running `install.packages("package_name")` for each package listed in `dependencies.R`.
5. Launch the application by running `shiny::runApp()`.

## Usage
0. Follow the instruction from the authore of qPCR analysis function:
 - The first column of the table should be named "sample" with names of samples and technical replicates.
 - The second column of the table should be named "group" with sample type. The rest cols are the CT values for each analyzed gene.
 - Import LinReg values with calculated primer amplification efficiency in xlsx.
 - To test the app use sample_data.xlsx
1. Choose the control group (for normalization):
2. Choose the vector with names of the reference genes (at least 2):
3. Load your PCR data using the provided import functionality.
4. Analyze the data as needed.
5. Explore the analysis options to uncover insights from your PCR results.
6. Visualize the data using various plotting functions.
7. Export your analysis results and visualizations for further use or sharing.

## Feedback
It is my first app using Shiny, created during my bachelor studies, so I am more than welcome to any feedback about the PCR Analyzer & Visualizer!

## Acknowledgements
This application was developed under the supervision of Ph. D. Piotr Gawroński and with the utilization of qPCR result analysis functions authored by Kinga Gołębiewska.
