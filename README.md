# EBBA_NichevsRangeShift_2022

Project investigating the niche and rage shifts represented between versions 1 and 2 of the European Breeding Bird Atlas.

##Current status
preliminary development of workflow and methods with raw EBBA1 and EBBA2 data. Many stages will be obsolete if and when the aligned EBBA1 and 2 data are available.

##Folder structure:
###scripts
Contains R scripts for; 
1. generating model inputs from raw EBBA files and environmental layers (Env_info_JH.md, Comparible_country_mask_generation.R, ecospat_input_preparation.R and species_list_generation.R)
2. Running niche and range shift analyses (ecospat_all_climate-change_run.R and ecospat_analysis_loop.R)
3. Combining species and model outputs into dataframes (ecospat_output_summation.R)
4. Exploratory plotting of combined species outputs (ecospat_combined_species_plots.R).

###plots 
Contains plots of finalized outputs. currently just for overall environmental similarity analysis.

###data (network archive only)
Contains raw data, and intermediary data and model outputs produced by scrips in the scripts folder.