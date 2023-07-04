# match phylogenetic information with our species names

library(dplyr)
library(phylolm)

# project data:
data_dir <- file.path("data")
data_dir_BBS <- file.path("data", "BBS_analysis")
data_dir_EBBA <- file.path("data", "EBBA_analysis")
plots_dir <- file.path("plots")

#-----------------------------------------------

# load traits
traits_all <- read.csv(file.path(data_dir,"Traits_niche_range_metrics.csv"))

# load phylogeny
load(file.path(data_dir,"PhyloBirdsHackett1_63.rdata"))

# identify non-matching species names
(no_match <- traits_all$species[!sub(' ','_',traits_all$species) %in% PhyloBirds$tip.label])

# add new columns 
traits_all$phylo_species <- traits_all$species
traits_all$phylo_species[traits_all$phylo_species=="Antrostomus carolinensis"] <- "Caprimulgus carolinensis"
traits_all$phylo_species[traits_all$phylo_species=="Antrostomus vociferus"] <- "Caprimulgus vociferus"
traits_all$phylo_species[traits_all$phylo_species=="Ardea alba"] <- "Casmerodius albus"
traits_all$phylo_species[traits_all$phylo_species=="Cardellina canadensis"] <- "Wilsonia canadensis"
traits_all$phylo_species[traits_all$phylo_species=="Cardellina pusilla"] <- "Wilsonia pusilla"
traits_all$phylo_species[traits_all$phylo_species=="Centronyx henslowii"] <- "Ammodramus henslowii"
traits_all$phylo_species[traits_all$phylo_species=="Circus hudsonius"] <- "Circus cyaneus"
# Careful: here two subspecies use the same tip - remove?
traits_all$phylo_species[traits_all$phylo_species=="Colaptes auratus auratus"] <- "Colaptes auratus"
traits_all$phylo_species[traits_all$phylo_species=="Colaptes auratus cafer"] <- "Colaptes auratus"
traits_all$phylo_species[traits_all$phylo_species=="Dryobates pubescens"] <- "Picoides pubescens"
traits_all$phylo_species[traits_all$phylo_species=="Dryobates scalaris"] <- "Picoides scalaris"
traits_all$phylo_species[traits_all$phylo_species=="Dryobates villosus"] <- "Picoides villosus"
traits_all$phylo_species[traits_all$phylo_species=="Gallinago delicata"] <- "Gallinago gallinago"
traits_all$phylo_species[traits_all$phylo_species=="Geothlypis formosa"] <- "Oporornis formosus"
traits_all$phylo_species[traits_all$phylo_species=="Geothlypis philadelphia"] <- "Oporornis philadelphia"
traits_all$phylo_species[traits_all$phylo_species=="Geothlypis tolmiei"] <- "Oporornis tolmiei"
traits_all$phylo_species[traits_all$phylo_species=="Haemorhous cassinii"] <- "Carpodacus cassinii"
traits_all$phylo_species[traits_all$phylo_species=="Haemorhous mexicanus"] <- "Carpodacus mexicanus"
traits_all$phylo_species[traits_all$phylo_species=="Haemorhous purpureus"] <- "Carpodacus purpureus"
# Careful: subspecies!
traits_all$phylo_species[traits_all$phylo_species=="Junco hyemalis hyemalis"] <- "Junco hyemalis"
traits_all$phylo_species[traits_all$phylo_species=="Junco hyemalis oreganus"] <- "Junco hyemalis"
traits_all$phylo_species[traits_all$phylo_species=="Leiothlypis celata"] <- "Vermivora celata"
traits_all$phylo_species[traits_all$phylo_species=="Leiothlypis ruficapilla"] <- "Vermivora ruficapilla"
traits_all$phylo_species[traits_all$phylo_species=="Mareca strepera"] <- "Anas strepera"
traits_all$phylo_species[traits_all$phylo_species=="Parkesia motacilla"] <- "Seiurus motacilla"
traits_all$phylo_species[traits_all$phylo_species=="Parkesia noveboracensis"] <- "Seiurus_noveboracensis"
traits_all$phylo_species[traits_all$phylo_species=="Peucaea cassinii"] <- "Aimophila cassinii"
traits_all$phylo_species[traits_all$phylo_species=="Poecile atricapillus"] <- "Parus atricapillus"
traits_all$phylo_species[traits_all$phylo_species=="Poecile carolinensis"] <- "Parus carolinensis"
traits_all$phylo_species[traits_all$phylo_species=="Poecile gambeli"] <- "Parus gambeli"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga americana"] <- "Parula americana"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga caerulescens"] <- "Dendroica caerulescens"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga cerulea"] <- "Dendroica cerulea"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga citrina"] <- "Wilsonia citrina"
# Careful: subspecies
traits_all$phylo_species[traits_all$phylo_species=="Setophaga coronata audoboni"] <- "Dendroica coronata"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga coronata coronata"] <- "Dendroica coronata"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga discolor"] <- "Dendroica discolor"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga dominica"] <- "Dendroica dominica"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga fusca"] <- "Dendroica fusca"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga magnolia"] <- "Dendroica magnolia"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga nigrescens"] <- "Dendroica nigrescens"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga pensylvanica"] <- "Dendroica pensylvanica"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga petechia"] <- "Dendroica petechia"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga pinus"] <- "Dendroica pinus"
traits_all$phylo_species[traits_all$phylo_species=="Setophaga virens"] <- "Dendroica virens"
traits_all$phylo_species[traits_all$phylo_species=="Spatula clypeata"] <- "Anas clypeata"
traits_all$phylo_species[traits_all$phylo_species=="Spatula discors"] <- "Anas discors"
traits_all$phylo_species[traits_all$phylo_species=="Spinus pinus"] <- "Carduelis pinus"
traits_all$phylo_species[traits_all$phylo_species=="Spinus psaltria"] <- "Carduelis psaltria"
traits_all$phylo_species[traits_all$phylo_species=="Spinus tristis"] <- "Carduelis tristis"
traits_all$phylo_species[traits_all$phylo_species=="Tringa semipalmata"] <- "Catoptrophorus semipalmatus"
traits_all$phylo_species[traits_all$phylo_species=="Troglodytes hiemalis"] <- "Troglodytes troglodytes"
traits_all$phylo_species[traits_all$phylo_species=="Vermivora cyanoptera"] <- "Vermivora chrysoptera"
traits_all$phylo_species[traits_all$phylo_species=="Emberiza calandra"] <- "Miliaria calandra"
traits_all$phylo_species[traits_all$phylo_species=="Linaria cannabina"] <- "Carduelis cannabina"
traits_all$phylo_species[traits_all$phylo_species=="Spinus spinus"] <- "Carduelis spinus"
traits_all$phylo_species[traits_all$phylo_species=="Cyanistes caeruleus"] <- "Parus caeruleus"
traits_all$phylo_species[traits_all$phylo_species=="Lophophanes cristatus"] <- "Parus cristatus"
traits_all$phylo_species[traits_all$phylo_species=="Periparus ater"] <- "Parus ater"
traits_all$phylo_species[traits_all$phylo_species=="Poecile lugubris"] <- "Parus lugubris"
traits_all$phylo_species[traits_all$phylo_species=="Poecile palustris"] <- "Parus palustris"
traits_all$phylo_species[traits_all$phylo_species=="Passer italiae"] <- "Passer hispaniolensis"

phylo_all= drop.tip(PhyloBirds,PhyloBirds$tip.label[!PhyloBirds$tip.label %in% sub(' ','_',traits_all$phylo_species)])
save(phylo_all,file=file.path(data_dir,"phylo_BBS_EBBA.RData"))

write.csv(traits_all, file=file.path(data_dir,"Traits_niche_range_metrics.csv"), row.names=F)


