# Range and climatic niche shifts in European and North American breeding birds
by Damaris Zurell, Katrin Schifferle, Sergi Herrando, Verena Keller, Aleksi Lehikoinen, Thomas Sattler, Levin Wiedenroth

**Summary**: Species respond dynamically to climate change leading to transient dynamics and time lagged range shifts. In consequence, species may not occupy their full climatic niche during range shifting. Here, we aim to assess climate niche tracking during recent range shifts of European and US birds. Using data from two European bird atlases and from the North American Breeding Bird Survey between 1980s–2010s, we analysed range overlap and climate niche overlap based on kernel density estimation and null model analyses. Phylogenetic multiple regression was used to assess the effect of species morphological, ecological and biogeographic traits on range and niche metrics. European birds shifted their ranges north and north-eastwards, US birds westwards. Range unfilling was much lower than expected by null models in both regions, suggesting inertia. Although transient dynamics were more pronounced in US, we found niche expansion was more common than niche unfilling in both regions. Trait analyses revealed some commonalities but also important differences between the two regions. Overall, dispersal limitations were minor in both regions while severe lagging of current range edges could be indicative of extinction debts. Differences between regions suggest diverging underlying mechanisms, for example related to land use history, that should be further investigated.

**Keywords**: range tracking, niche tracking, transient dynamics, extinction debts, colonisation credit, dispersal limitation

This repository contains the R scripts needed to reproduce all results as well as results and plots. Codes were implemented by Katrin Schifferle with contributions from Damaris Zurell.

**Funding**: This research was supported by Deutsche Forschungsgemeinschaft (DFG) under grant agreements No. ZU 361/1-1 and ZU 361/6-1.

---

## General workflow:

### EBBA analysis:

In a first step, EBBA1 and EBBA2 data were processed into spatial data (1_EBBA_prep_data.R). Then a subset of the species considered suitable for the analyses was extracted: pelagic specialists and very rare and common species were excluded (2_1_EBBA_species_filtering_1-4.R, final_filtering = FALSE), afterwards species whose climatic niche is not well covered in Europe (stability < 50 %) were excluded (2_3_EBBA_species_filtering_5_climatic_niche_analysis.R). Left were 118 species for which the EBBA1-EBBA2 change dataset was requested. This dataset includes only the EBBA cells for which the data are comparable across both atlas versions. The EBBA1-EBBA2 change dataset was then converted into spatial data (1_EBBA_prep_data.R). These data were again filtered to exclude very rare and common species now based on the EBBA1-EBBA2 change dataset (2_1_EBBA_species_filtering_1-4.R, final_filtering = TRUE). Left were 117 species for which then the range and niche shift analyses were conducted (4_EBBA_range_shift_analysis.R, 4_EBBA_niche_shift_analysis.R). To quantify the climatic niche, bioclimatic variables were calculated for Europe (3_EBBA_prep_env_data.R). Results of the range and niche shift analyses were explored and plotted (5_EBBA_plots_results.R). For a pending analysis on whether species traits influence niche and range dynamics, species traits based on AVONET were assembled and the global niche breadths were quantified (3_prep_trait_data.R). Exploratory plots regarding traits and niche and range dynamics were made (5_trait_plots.R).

### BBS analysis:
In a first step, BBS data were processed into spatial presence/absence data (1_BBS_prep_data.R). Only the BBS routes that were sampled in all three years of both the historic and the recent time period were used.
Then a subset of the species considered suitable for the analyses was extracted: pelagic specialists, very rare and common species as well as records for which data were not available on the species level were excluded (2_1_BBS_species_filtering_1-4.R), afterwards species whose climatic niche is not well covered in the conterminous US (stability < 50 %) were excluded (2_3_BBS_species_filtering_5_climatic_niche_analysis.R). Left were 207 species when using 1981-1983 as the historic time period and 296 species when using 1996-1998 as the historic time period. For these species the range and niche shift analyses were conducted (4_BBS_range_shift_analysis.R, 4_BBS_niche_shift_analysis.R). To quantify the climatic niche, bioclimatic variables that were calculated for the conterminous US (3_BBS_prep_env_data.R). Results of the range and niche shift analyses were explored and plotted (5_BBS_plots_results.R). For a pending analysis on whether species traits influence niche and range dynamics, species traits based on AVONET were assembled and the global niche breadths were calculated (3_prep_trait_data.R). Exploratory plots regarding traits and niche and range dynamics were made (5_trait_plots.R).


## Abbreviations in file names:

- EBBA = analysis for European breeding bird species based on the European Breeding Bird Atlas (EBBA)
- BBS = analysis for North American breeding bird species based on the North American Breeding Bird Survey (BBS)
- BL = BirdLife 2022 range maps

### bg_spec / bg_EBBA / bg_US:

Niche and range shift analysis were run using two versions of environmental backgrounds, i.e. the total area that is assumed to be available to a species:
- "bg_spec" = species-specific environmental background: for Europe we use the EBBA-change cells located within a 500 km buffer around the species presences, for the conterminous US we use the routes with presences and true absences of a species within a 500 km buffer around the presences (following Sofaer et al.)
- "bg_EBBA" / "bg_US" = same environmental background for all species: for Europe we use all EBBA-change cells, for the US we use all routes containing presences and true absences of a species across the conterminous US

### hist81-83 / hist96-98 (BBS only):

The BBS analysis was run separately for two different historic time periods:

- "hist81-83" = 1981 - 1983: earliest time span for which Chelsa data is available / maximum gap between historic and recent time period (Chelsa data available from 1980 onwards, we included the year preceding the considered time period in calculating the bioclimatic variables to account for lag effects following Sofaer et al.)
- "hist88-90" = 1988 - 1990: yields a similar time gap between the historic and the recent time period as in the EBBA analysis


## Folder structure:

### data:

- input data necessary to run the scripts that is currently not stored on the labs' datashare (EltonTraits: BirdFuncDat.txt, AVONET trait data: AVONET Supplementary dataset 1.xlsx)

#### BBS_analysis:

- some of the smaller output files produced by running the BBS analysis scripts

#### EBBA_analysis:

- some of the smaller output files produced by running the EBBA analysis scripts


### plots: 

- niche_breadth_species: climatic niche space occupied by the species compared to global climatic niche space, plot titles: species and niche breadth (output of: 3_prep_trait_data.R)
- coverage_climatic_niche: plots regarding niche coverage in Europe / in the conterminous US: comparison of climatic niche used throughout the global range based on Birdlife 2022 maps (red lines) and climatic niche covered by the EBBA 2 occurrences / by the part of the global range located in the conterminous US (green lines; underlying density corresponds to global range; plots are output of 2_3_species_filtering_5_climatic_niche_analysis.R)
- niche_dynamics_species: niche dynamics (stability = blue, expansion = red, unfilling = green) of each species in 2-dimensional niche space (output of 4_niche_shift_analysis.R)
- range_dynamics_species: range dynamics (stability = blue, expansion = red, unfilling = green) of each species (output of 4_range_shift_analysis.R)
- dynamics_correlations: scatterplots of values for range and niche dynamics to explore correlations (output of 5_plots_results.R)
- dynamics_boxplots: boxplots of values for range and niche dynamics (output of 5_plots_results.R)
- climate_change_PCAs: explore climate change in Europe / the conterminous US between the historic and the recent time period: PCA based biplot for bioclimatic variables of EBBA area / conterminous US and maps showing the differences in PC1 and PC2 between both time periods (output of 5_plots_results.R)
- range_shift_direction: directions of the shift between the range centroid of the historic period and the range centroid of the recent period (output of 5_plots_results.R)
- species_richness: maps of the distribution of the number of species recorded in the comparable EBBA cells / the selected conterminous US routes in both time periods (output of 5_plots_results.R)
- traits_explorations: exploratory plots regarding the relationship between species traits and range / niche dynamics (overlap D, expansion E, stability S, unfilling U) (output of 5_trait_plots.R)


### scripts:

Script file names that don't contain "EBBA" or "BBS" work for both analyses. In the scripts' set-up section it can be specified for which analysis (EBBA or BBS) the script should be run.
It can also be specified which version of environmental background and which historic time period (for the BBS analysis) should be used and file paths can be adjusted.

#### 1_(EBBA/BBS)_prep_data.R:

##### EBBA:
Preparations of EBBA 1 and 2 data: 
- projecting spatial features using Lambert azimuthal equal-area projection (ETRS89-extended / LAEA Europe, EPSG:3035)
- dropping EBBA cells that are not comparable across both atlas versions (outputs: EBBA1_prelim_comparable_cells.shp and EBBA2_prelim_comparable_cells.shp)
- dealing with species name changes between EBBA 1 and EBBA 2 (outputs: EBBA1_prelim_comparable_harmonized.shp, EBBA2_prelim_comparable_harmonized.shp)

Processing of EBBA1-EBBA2 change dataset (outputs: EBBA_change.shp, EBBA1_change.shp, EBBA2_change.shp)

##### BBS:
- filtering BBS routes: keep routes that meet official BBS criteria and that were sampled in all years of both time periods, remove water based routes
- converting counts to presence/true absence data (presence = a species was present at least once in a time period on a route, true absence = a species was not observed in any year of the time period, although the route was sampled in all three years)
- calculating route centroids if full spatial information of route is available, use route starting point if not
- projecting spatial features using Albers equal-area projection (ESRI:102003)
- outputs: BBS_historic_centr_proj_hist81-83.shp, BBS_historic_centr_proj_hist96-98.shp, BBS_recent_centr_proj_hist81-83.shp, BBS_recent_centr_proj_hist96-98.shp

#### 2_1_(EBBA/BBS)_species_filtering_1-4.R:

##### EBBA:
Subset the species to keep only those that are probably most suitable for the analysis:
- step 1: Exclude pelagic specialists (according to Wilman et al. 2014)
- step 2: Exclude rare species with n<20 occurrences in any of the two atlas periods.
- step 3: Exclude very common species with >90% prevalence in at least one of the two atlas periods.
- step 4: Keep only those species that are native occur in both atlas periods.
- outputs: EBBA1_EBBA2_prep_steps1-4_prelim.RData, EBBA1_EBBA2_prep_steps1-4_final.RData

##### BBS:
Same as for EBBA, but also excluding records for which data were not available on the species level.
- outputs: BBS_prep_steps1-4_hist81-83.RData, BBS_prep_steps1-4_hist88-90.RData

#### 2_2_species_filtering_5_project_Chelsa.R:

Preparation for step 5 to further subset the species: Chelsa data of 2012 - 2018 (including the time period covered by EBBA 2 and the recent period of the BBS analysis) are projected using the global equal area projection <em>Interrupted Goode Homolosine</em> and a resolution of 50 km, which matches the resolution of the EBBA.

#### 2_3_(EBBA/BBS)_species_filtering_5_climatic_niche_analysis.R:

##### EBBA:

To exclude species whose climatic niche is not well covered in Europe, we compared the climatic niche that a species uses throughout its global range with the climatic niche covered by the comparable EBBA area. To characterise the climatic niche we use 19 bioclimatic variables derived from Chelsa data of the years 2012 - 2017.
- to describe the climatic niche, occurrence records are matched to the corresponding bioclimatic variables: with regard to the European distribution the EBBA 2 occurrence records are used, with regard to the global distributions the raster cell centroids of the Birdlife 2022 range that is either used throughout the year or during the breeding season are used as occurrence points
- following the niche quantification and comparison workflow of the <em>ecospat</em> package, the climatic niche of each species is summarised using the first two axes of a PCA. Afterwards occurrence density grids are calculated in climatic niche space for both niches (global as z2 and EBBA-based as z1) and niche dynamics are quantified. The stability index is used to determine how much of the species' global climatic niche is covered in Europe. 
- in the following analyses we used the EBBA1-EBBA2 change dataset for species with stability >= 50 % (previously, we used the 150 species with the largest stability values, thus requested in fact EBBA1-EBBA2 change data for some more species)
- outputs: species_stability_EBBA2_BL22_030423.csv, Bioclim_global_2009_2018, plots/coverage_climatic_niche/EBBA2_vs_BL_niche_dyn

##### BBS:

Similar as above, however, instead of comparing the BBS distribution of a species to its global distribution, we assessed how much of a species global climatic niche is contained in the conterminous US:
- Chelsa data of the years 2015 - 2018 are used
- as occurrences records we used the raster cell centroids of the Birdlife 2022 range that is either used throughout the year or during the breeding season
- niche quantification: global climatic niche as z2 and climatic niche contained in the conterminous US as z1
- outputs: species_stability_contUS_BL22.csv, Bioclim_global_2015_2018, plots/coverage_climatic_niche/contUS_vs_BL_niche_dyn, conterminousUS.shp

#### 3_(EBBA/BBS)_prep_env_data.R:

##### EBBA:
Chelsa climate data of the area covered by the EBBA 2 grid (mask: EBBA2_area.tif) of the historic (1984-1988) and recent (2012-2017) time period are projected using the Lambert azimuthal equal-area projection (ETRS89-extended / LAEA Europe, EPSG:3035) and a resolution of 50km. Subsequently, the values of the 19 bioclimatic variables are calculated based on the monthly mean values of precipitation, minimum, maximum and mean temperature, once for the historic and once for the recent time period.
- outputs: Bioclim_1984_1988, Bioclim_2012_2017

##### BBS:
Chelsa climate data of the conterminous US (mask: conterminousUS_1km.tif) of the historic and recent time periods including the preceding year to account for lag effect (see Sofaer et al.; 1980-1883 / 1987-1990) are projected using the Albers equal-area projection (ESRI:102003) and a resolution of 1 km. Subsequently, the values of the 19 bioclimatic variables are calculated as above.
- outputs: Bioclim_1980_1983, Bioclim_1987_1990, Bioclim_2015_2018


#### 3_prep_trait_data.R:
- species that passed the filtering steps are matched with AVONET trait database
- the global climatic niche breadth of each species is quantified by first determining the climatic niche space of the species global range with the <em>ecospat</em> package, using the same, global environmental background for all species to get comparable results, and then calculating the Shannon index of the resulting occurrence density grid corrected for environmental prevalence as a measure of the global niche breadth
- outputs: AVONET_(EBBA/BBS)_species.csv, (EBBA/BBS)_niche_breadth.csv

#### 4_(EBBA/BBS)_niche_shift_analysis.R:

Niche overlap analysis using the <em>ecospat</em> package. The niche is calculated as the density scores along the first two axes of a PCA that summarises the bioclimatic variables. 
(For the EBBA analysis the values of the bioclimatic variables at the cell centroids are used. For the BBS analysis the climate data corresponding to a route centroid are determined by first buffering the routes' centroids with a 21 km buffer and then calculting a weighted mean of the values of the bioclimatic variables within the buffer, see Sofaer et al.). For each species Schoener's D, the results of niche similarity tests with regard to niche shifting (column names starting with "shift_p") and niche conservatism (column names starting with "cons_p") with regard to analogue (column name containing "A") and non-analogue (column name containing "NA") conditions, as well as standardized indices for niche dynamics (stability S, expansion E, unfilling U, abandonment A and pioneering P) are stored in a csv file.
- outputs: EBBA_niche_shift_results_bg_EBBA.csv, EBBA_niche_shift_results_bg_spec.csv, BBS_niche_shift_results_bg_spec_hist81-83.csv, BBS_niche_shift_results_bg_US_hist81-83.csv, BBS_niche_shift_results_bg_spec_hist88-90.csv, BBS_niche_shift_results_bg_US_hist88-90.csv, plots/niche_dynamics_species


#### 4_(EBBA/BBS)_range_shift_analysis.R:

First, for each species the central range coordinates in the historic and recent time period (columns "centroid_hist_X", "centroid_hist_Y", "centroid_rec_X", "centroid_rec_Y"), the Euclidian distance (column "Eucl_distance") between them, the shift in north-south (column "NS_shift") and east-west direction (column "EW_shift"), the overall direction of the shift (column "shift_direction") as well as the respective range size (columns "n_cells_hist", "n_cells_rec") is calculated. Secondly, a range overlap analysis is conducted analogous to the niche analysis using the <em>ecospat</em> package. Instead of conducting a PCA to reduce the environmental data to 2 dimensions, the analysis is conducted in geographic space using the X and Y coordinates rather than the first two PCA axes. 
The density grids of species occurrences are restricted to the EBBA land area (EBBA_ecospat_mask_change_cells_buffer20km.tif) and the conterminous US, Canada and Mexico, respectively. Schoener's D, the results of similarity tests with regard to range shifting (column names starting with "shift_p") and range conservatism (column names starting with "cons_p") with regard to analogue (column name containing "A") and non-analogue (column name containing "NA") conditions, as well as standardized indices for range dynamics (stability S, expansion E and unfilling U) are stored in a csv file.
- outputs: EBBA_range_shift_results_bg_EBBA.csv, EBBA_range_shift_results_bg_spec.csv, BBS_range_shift_results_bg_spec_hist81-83.csv, BBS_range_shift_results_bg_US_hist81-83.csv, BBS_range_shift_results_bg_spec_hist96-98.csv, BBS_range_shift_results_bg_US_hist96-98.csv, plots/range_dynamics_species

#### 4_trait_analysis.R and 4_prep_phylo_data.R
Scripts to prepare the trait and phylogenetic data and run phylogenetic regression models, incl. permutation approach to quantify variable importance.

#### 5_(EBBA/BBS)_plots_results.R and 5_joint_plots_results.R:
Scripts to produce the plots (see "plots" section).

#### 5_trait_plots.R:
Script to produce exploratory plots regarding the relationship between species traits and range / niche dynamics (see "plots" section).
