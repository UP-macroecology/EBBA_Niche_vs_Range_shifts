Scripts;

1. env_info_JH.Rmd - Adapted version of script produced by Simon Kapitza and Katrin Schifferle. Produces Bioclimatic and Land-Use layers for project at suitable extent and resolution (from CHELSA and LUH2 sources).

2. Comparable country mask.R - Produces country masks from an international boundary shapefile, suitable for rough spatial limitation of EBBA and climate data to a comparable region indicated in EBBA_methods chapter. Output is two shapefiles, MaskCountries and KeepCountries, used in following scripts.

3. ecospat_input_preparation.R - Restricts Bioclimatic data for ecospat analyses (masking to suitable study region).  

4. ecospat_all_climate_change_run.R - Full environment run of ecospat analyses to assess overall level of climate change between study periods. Outputs: Shoeners'D and other metrics, niche similarity test results, kernel density plot and niche similarity test plot.

5. Species_list_generation.R - Taxonomic unification between EBBA1 and EBBA2 data. Derivation of species list for request for aligned dataset: steps: remove species with <20 occurrences in either EBBA1 or EBBA2. Plot all distributions to assess any remaining issues visually. Remove pelagic species based on EltonTraits database. Remove species recorded in >90% comparable EBBA cells. Remove species with less than 5% overlap between IUCN range polygon and comparable EBBA cells. The final step is likely not defensible or sensible and needs reconsidering. Outputs: Provisional species list for submission for aligned EBBA data. Filtered EBBA1 and EBBA2 datasets to be used in preliminary model runs.

6. ecospat_analysis_loop.R - Calculates niche and range shifts from EBBA data using ecospat package. Main analyses conducted in foreach loop with outputs written to file. Outputs, for each species: PCA axes contributions, comparison of presence and absence cells in EBBA 1 and 2, Range size calculation for each EBBA version, range expansion, contraction and stability metrics, range centroids and range shift metric, output from ecospat analysis for 5 buffer sizes and 4 versions of niche similarity test, kernel density and niche similarity plots.

7. ecospat_output_summation.R - Coerces output files from ecospat_analysis_loop.R into summary dataframes that can be analysed and plotted further. Outputs include several dataframes stored as .csv files. Included relative dynamic calculations provided by AR.


8. ecospat_combined_species_plots.R - Produces exploratory plots of combined niche and range shift outputs. Outputs; buffer size sensitivity analysis, niche and range shift correlations (grouped and ungrouped) and expansion, contraction and stability metric correlations in real and niche space.
