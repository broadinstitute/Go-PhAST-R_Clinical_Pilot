# Go-PhAST-R_Clinical_Pilot

Scripts:
- HeatMapFuncs.R
  - Accessory functions for NSTGdata_calculations.R
- NSTGdata_calculations.R - Normalization and calculation script
  - Normalization to positive and negative controls
  - Normalization to control genes
  - Absolute average probe response data calculations
  - Log2 fold change between treated and untreated calculations
  - Squared projected distance calculations
- Go-PhAST-R_vis.Rmd - Visualization and prediction script
  - Create figures
  - SVM

Files:
- Raw data .csv files
- Log2 fold change .csv files
- Normalized beta lactamase .csv file
- mel data here

Use:

NSTGdata_calculations.R
- This script expects raw data in the same format as the raw data .csv files:
    1. Columns - by sample, named, organized by clinical strain, alternating untreated first and then treated
    2. Rows - probe response, named
    3. giving a directory name with the files [dfnames] is sufficient to import each of these
- [known_data_in] - This script also requires a dataset from strains with known MICs. This data is 'long' form meaning that there are columns for each Strain, drug treatment, probe response value, SIR, and MIC.
- [CRE_bckgrnd] - This file has the background signal for each beta lactamase in our assay based on the average signal plus two standard deviations of previous data confirmed to be negative for that geneotype.
Adjust paths and file names for input and output and run script.

Go-PhAST-R_vis.Rmd
- This




