# Patterns of niche and range shifts in European breeding birds under climate change

Project investigating niche and range shifts in European breeding birds by comparing versions 1 and 2 of the European Breeding Bird Atlas (EBBA).

## Current status

The EBBA 1 and EBBA 2 cells that are comparable across both atlas versions according to the EBBA 2 methods chapter are stored as shapefiles.
Then a subset of these shapefiles was created containing only the species for which the two atlas versions can be compared and in which the species names of EBBA 1 were adjusted to the names used in EBBA 2.
A subset of 150 species was determined that is probably most suited for the analysis. For these species we request the EBBA1-EBBA2 change dataset.

## Folder structure:

### Data:

- input data necessary to run the scripts that is currently not stored on the labs' datashare (BirdFuncDat.txt, ranges of <em>Passer italiae</em> and <em>Sylvia ruppeli</em>)
- some of the smaller output files produced by running the scripts (Bioclim_raster_EBBA1, Birdlife_ranges_mask.tif, EBBA1_EBBA2_prep_steps1-4.RData, EBBA_niche_range_shifts_species_selection.csv)

### plots: 

- EBBA-Birdlife_niche_dyn: plots regarding niche coverage in Europe for each of the 325 species left after filtering step 4 (see below): comparison of climatic niche used throughout the global range based on Birdlife maps (red lines) and climatic niche covered in Europe (green lines; underlying density corresponds to global range; plots are output of 2_3_species_filtering_5_climatic_niche_analysis.R)

### scripts:

#### 1_prep_EBBA_data.R: 

Preparations of EBBA 1 and 2 data: 
- projecting spatial features using Lambert azimuthal equal-area projection (ETRS89-extended / LAEA Europe, EPSG:3035)
- dropping EBBA cells that are not comparable across both atlas versions (outputs: EBBA1_comparable.shp and EBBA2_comparable.shp)
- dealing with species name changes between EBBA 1 and EBBA 2: if species records can be compared, the name that is used in EBBA 2 is also used for EBBA 1 records, if a taxa was treated as a single species in EBBA 1 but split into several species in EBBA 2, only the species that can be compared nevertheless according to the EBBA 2 methods chapter are retained (outputs: "EBBA1_comparable_harmonized.shp", "EBBA2_comparable_harmonized.shp")

#### 2_1_species_filtering_1-4.R:

Steps 1 to 4 to subset the species to keep only those that are probably most suitable for the analysis:
- inputs: EBBA1_comparable_harmonized.shp, EBBA2_comparable_harmonized.shp, Data/BirdFuncDat.txt
- step 1: Exclude pelagic specialists (according to Wilman et al. 2014)
- step 2: Exclude rare species with n<20 occurrences in any of the two atlas periods.
- step 3: Exclude very common species with >90% prevalence in at least one of the two atlas periods.
- step 4: Keep only those species that occur in both atlas periods.
- output: Data/EBBA1_EBBA2_prep_steps1-4.RData

#### 2_2_species_filtering_5_project_Chelsa.R:

Preparation for step 5 to further subset the species: Chelsa data of the time period covered by EBBA 1 (1981 - 1990) are projected using the global equal area projection <em>Interrupted Goode Homolosine</em> and a resolution of 50 km, which matches the resolution of the EBBA.

#### 2_3_species_filtering_5_climatic_niche_analysis.R:

Step 5 to subset the species: exclude species whose climatic niche is not well covered in Europe. 
We compared the climatic niche that a species uses throughout its global range with the climatic niche covered by the comparable EBBA area. To characterise the climatic niche we use 19 bioclim variables derived from Chelsa data.
- inputs: Data/EBBA1_EBBA2_prep_steps1-4.RData, EBBA1_comparable.shp
- A mask is created to extract only the relevant Chelsa data, it covers the global ranges of all 325 species left after filtering step 4, the area where EBBA 1 and EBBA 2 are comparable and a buffer of 100 km. The mask is then applied to the Chelsa data and from the masked Chelsa data monthly means with regard to temperature and precipitation are determined, which in turn are used to calculate the 19 bioclim variabels (output: Data/Bioclim_raster_EBBA1/...).
- As a basis to describe climatic niches, occurrence records of each species are matched to the corresponding bioclim variables: with regard to the European distribution the EBBA 1 occurrence records are used, with regard to the global distributions the raster grid points of the Birdlife range that is either used throughout the year or during the breeding season are used as occurrence points.
- Following the niche quantification and comparison workflow of the <em>ecospat</em> package, a PCA of the resulting datasets is conducted and the first 2 PCA axes are used to describe the climatic niche of each species. Afterwards Occurrence Density Grids are calculated in climatic niche space for both niches (global as z2 and EBBA-based as z1) and niche dynamics are quantified. The stability index is used to quantify how much of the species' global climatic niche is covered in Europe. 
- For now, we selected the 150 species with the largest stability values. For these species the EBBA1-EBBA2 change dataset is requested.
- outputs: Data/EBBA_niche_range_shifts_species_selection.csv, plots/EBBA-Birdlife_niche_dyn/...

The folder <em>deprecated</em> contains a preliminary workflow for the whole project.