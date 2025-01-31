# Go-PhAST-R_Clinical_Pilot

Young, E.L., Roach, D.J., Martinsen, M.A., McGrath, G.E., Holbrook, N.R., Cho, H.E., Seyoum, E.Y., Pierce, V.M. and Bhattacharyya, R.P., 2024. Clinical pilot of bacterial transcriptional profiling as a combined genotypic and phenotypic antimicrobial susceptibility test. Journal of Clinical Microbiology, 62(11), pp.e00997-24.

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
- 2021.12.13_CRE_bkgd_subt_mean2SD.csv: Background signal for absent beta lactamase genes calculated from samples with known genotypes.
- AllMelpaperData_all.csv: Data including log2 fold change, strain, MIC, and drug treatment from Martinsen MA, Jaramillo Cartagena A, Bhattacharyya RP. Core Antibiotic-Induced Transcriptional Signatures Reflect Susceptibility to All Members of an Antibiotic Class. Antimicrob Agents Chemother. 2021;65(6):10.1128/aac.02296-20. doi:10.1128/aac.02296-20
- ClinData.csv: This data includes the SIR call and MIC for each strain and drug combination as determined by the clinical microbiology lab.
- probe_to_gene_names.csv: This .csv matches NanoString probe names with the target gene.

Use:

NSTGdata_calculations.R
- This script expects raw data in the same format as the raw data .csv files:
    1. Columns - by sample, named, organized by clinical strain, alternating untreated first and then treated
    2. Rows - probe response, named
    3. giving a directory name with the files [dfnames] is sufficient to import each of these
- [known_data_in] - This script also requires a dataset from strains with known MICs. This data is 'long' form meaning that there are columns for each Strain, drug treatment, probe response value, SIR, and MIC.
- [CRE_bckgrnd] - This file has the background signal for each beta lactamase in our assay based on the average signal plus two standard deviations of previous data confirmed to be negative for that geneotype.
Adjust paths and file names for input and output and run script. It will generate several files incuding the following which are passed on to the visualization script:
- NSTGsample_BLase.csv: Binary response data for each strain and beta lactamase probe indicating presence and absence of the genotype.
- all_mel_data_SPD.csv: SPD data calculated from AllMelpaperData_all.csv .
- MGH_AST_spd.csv: SPD calculated for each sample.
- NSTGsample_CREnorm_processed.csv: Normalized count data for genotypic NanoString probe from each beta lactam sample.


Go-PhAST-R_vis.Rmd





