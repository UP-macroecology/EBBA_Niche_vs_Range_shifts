# Trait analyses using phylogenetic regression models
# Trait models are fitted separately to BBS and EBBA data as the niche breadth was calculated in different years?

library(dplyr)
library(phylolm)
# library(caper)


# project data:
data_dir <- file.path("data")
data_dir_BBS <- file.path("data", "BBS_analysis")
data_dir_EBBA <- file.path("data", "EBBA_analysis")
plots_dir <- file.path("plots")

#-------------------------------------------------------------------
# prepare trait data and range and niche metrics

# Niche and range analyses results for EBBA species
niche_results_EBBA <- read.csv(file.path(data_dir_EBBA, "EBBA_niche_shift_results_bg_spec_070723.csv")) %>%
  rename(niche_D = D) %>%
  select(species, niche_D, niche_stability_std, niche_unfilling_std, niche_expansion_std)
range_results_EBBA <- read.csv(file.path(data_dir_EBBA, "EBBA_range_shift_results_bg_spec_070723.csv")) %>%
  rename(range_D = D, range_shift = Eucl_distance) %>%
  select(species, range_shift, range_D, range_stability_std, range_unfilling_std, range_expansion_std)
# EBBA species list
sel_species_EBBA <- niche_results_EBBA$species

# Niche and range analyses results for EBBA species
niche_results_BBS <- read.csv(file.path(data_dir_BBS, "BBS_niche_shift_results_bg_spec_hist81-83_070723.csv")) %>%
  rename(niche_D = D) %>%
  select(species, niche_D, niche_stability_std, niche_unfilling_std, niche_expansion_std)
range_results_BBS <- read.csv(file.path(data_dir_BBS, "BBS_range_shift_results_bg_spec_hist81-83_070723.csv")) %>%
  rename(range_D = D, range_shift = Eucl_distance) %>%
  select(species, range_shift, range_D, range_stability_std, range_unfilling_std, range_expansion_std)
# BBS species list
sel_species_BBS <- niche_results_BBS$species

# We exclude species that occur in both regions (n=1)
remove_spp <- intersect(sel_species_EBBA, sel_species_BBS)
sel_species_EBBA <- sel_species_EBBA[!sel_species_EBBA %in% remove_spp]
sel_species_BBS <- sel_species_BBS[!sel_species_BBS %in% remove_spp]

# Trait data for EBBA species
avonet_EBBA <- read.csv(file.path(data_dir_EBBA, "AVONET_EBBA_species.csv")) %>%
  filter(Species1 %in% sel_species_EBBA)
niche_breadth_EBBA <- read.csv(file.path(data_dir_EBBA, "EBBA_niche_breadth_060723.csv")) %>%
  filter(species %in% sel_species_EBBA)

# Trait data for BBS species
avonet_BBS <- read.csv(file.path(data_dir_BBS, "AVONET_BBS_species.csv")) %>%
  filter(BBS_species %in% sel_species_BBS)
niche_breadth_BBS <- read.csv(file.path(data_dir_BBS, "BBS_niche_breadth_060723.csv")) %>%
  filter(species %in% sel_species_BBS)

# Combine trait data sets and select relevant traits
 traits_BBS <- inner_join(avonet_BBS, niche_breadth_BBS, join_by(BBS_species==species)) %>%
  rename(species = BBS_species) %>%
  select(c(species, Mass, Hand.Wing.Index, Trophic.Level, Habitat.Density, Range.Size, Migration, Centroid.Latitude, Centroid.Longitude, niche_breadth_zcor)) %>%
  mutate(region="US") %>%
   inner_join(range_results_BBS) %>%
   inner_join(niche_results_BBS)

traits_EBBA <- inner_join(avonet_EBBA, niche_breadth_EBBA, join_by(Species1==species)) %>%
  rename(species = Species1) %>%
  select(c(species, Mass, Hand.Wing.Index, Trophic.Level, Habitat.Density, Range.Size, Migration, Centroid.Latitude, Centroid.Longitude, niche_breadth_zcor)) %>%
  mutate(region="Europe") %>%
  inner_join(range_results_EBBA) %>%
  inner_join(niche_results_EBBA)

traits_all <- rbind(traits_BBS, traits_EBBA)
# Trophic level as ordinal scale from 1=herbivore to 3=carnivore & scavenger
traits_all$Trophic.Level <- ifelse(traits_all$Trophic.Level=="Carnivore",3,ifelse(traits_all$Trophic.Level=="Scavenger",3,ifelse(traits_all$Trophic.Level=="Omnivore",2,1)))

write.csv(traits_all, file=file.path(data_dir,"Traits_niche_range_metrics.csv"), row.names=F)

# switch to 4_prep_phylo_data.R to match names with phylogeny

#--------------------------------------------------

#load phylogenetic tree (named "phylo_all")
load(file.path(data_dir,"phylo_BBS_EBBA.RData"))

#--------------------------------------------------

# Read in trait data
traits_all <- read.csv(file.path(data_dir,"Traits_niche_range_metrics.csv"))

# remove subspecies (duplicate phylonames)
traits_all <- traits_all[!traits_all$phylo_species %in% traits_all$phylo_species[duplicated(traits_all$phylo_species)],]

# also remove from phylogeny
phylo_all= drop.tip(phylo_all,phylo_all$tip.label[!phylo_all$tip.label %in% sub(' ','_',traits_all$phylo_species)])

# add rownames - needed for matching phylogenetic information
rownames(traits_all) = sub(' ','_',traits_all$phylo_species)

# Define logit function
logit = function(x) {x=ifelse(x<0.0001,0.0001,ifelse(x>0.9999,.9999,x));log(x/(1 - x))}

# Standardise traits
traits_all$Mass <- as.numeric(scale(traits_all$Mass))
traits_all$Hand.Wing.Index <- as.numeric(scale(traits_all$Hand.Wing.Index))
traits_all$Range.Size <- as.numeric(scale(traits_all$Range.Size))
traits_all$Centroid.Latitude <- traits_all$Centroid.Latitude/90
traits_all$Centroid.Longitude <- traits_all$Centroid.Longitude/180
traits_all$niche_breadth_zcor <- as.numeric(scale(traits_all$niche_breadth_zcor))
traits_all$range_shift <- as.numeric(scale(traits_all$range_shift))
traits_all$region <- as.factor(traits_all$region)

# Logit-transform response variables
traits_all$range_D_logit <- logit(traits_all$range_D)
traits_all$range_stability_std_logit <- logit(traits_all$range_stability_std)
traits_all$range_unfilling_std_logit <- logit(traits_all$range_unfilling_std)
traits_all$range_expansion_std_logit <- logit(traits_all$range_expansion_std)
traits_all$niche_D_logit <- logit(traits_all$niche_D)
traits_all$niche_stability_std_logit <- logit(traits_all$niche_stability_std)
traits_all$niche_unfilling_std_logit <- logit(traits_all$niche_unfilling_std)
traits_all$niche_expansion_std_logit <- logit(traits_all$niche_expansion_std)


#-------------- MODELS ---------------

# range_D 
# null model
null.rangeD = phylolm(range_D_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.rangeD = phylolm(range_D_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor + region, 
                    data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.rangeD <- phylostep(lm.rangeD$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.rangeD$residuals^2)/sum(null.rangeD$residuals^2)
1-sum(step.lm.rangeD$residuals^2)/sum(null.rangeD$residuals^2)



# range_stability_std_logit 
# null model
null.range_stab = phylolm(range_stability_std_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.range_stab = phylolm(range_stability_std_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor + region, 
                    data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.range_stab <- phylostep(lm.range_stab$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.range_stab$residuals^2)/sum(null.range_stab$residuals^2)
1-sum(step.lm.range_stab$residuals^2)/sum(null.range_stab$residuals^2)



# range_unfilling_std_logit 
# null model
null.range_unf = phylolm(range_unfilling_std_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.range_unf = phylolm(range_unfilling_std_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor + region, 
                        data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.range_unf <- phylostep(lm.range_unf$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.range_unf$residuals^2)/sum(null.range_unf$residuals^2)
1-sum(step.lm.range_unf$residuals^2)/sum(null.range_unf$residuals^2)


# range_expansion_std_logit 
# null model
null.range_exp = phylolm(range_expansion_std_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.range_exp = phylolm(range_expansion_std_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor + region, 
                       data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.range_exp <- phylostep(lm.range_exp$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.range_exp$residuals^2)/sum(null.range_exp$residuals^2)
1-sum(step.lm.range_exp$residuals^2)/sum(null.range_exp$residuals^2)


# niche_D 
# null model
null.nicheD = phylolm(niche_D_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.nicheD = phylolm(niche_D_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor + region, 
                    data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.nicheD <- phylostep(lm.nicheD$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.nicheD$residuals^2)/sum(null.nicheD$residuals^2)
1-sum(step.lm.nicheD$residuals^2)/sum(null.nicheD$residuals^2)



# niche_stability_std_logit 
# null model
null.niche_stab = phylolm(niche_stability_std_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.niche_stab = phylolm(niche_stability_std_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor + region, 
                    data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.niche_stab <- phylostep(lm.niche_stab$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.niche_stab$residuals^2)/sum(null.niche_stab$residuals^2)
1-sum(step.lm.niche_stab$residuals^2)/sum(null.niche_stab$residuals^2)



# niche_unfilling_std_logit 
# null model
null.niche_unf = phylolm(niche_unfilling_std_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.niche_unf = phylolm(niche_unfilling_std_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor + region, 
                        data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.niche_unf <- phylostep(lm.niche_unf$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.niche_unf$residuals^2)/sum(null.niche_unf$residuals^2)
1-sum(step.lm.niche_unf$residuals^2)/sum(null.niche_unf$residuals^2)



# niche_expansion_std_logit 
# null model
null.niche_exp = phylolm(niche_expansion_std_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.niche_exp = phylolm(niche_expansion_std_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor + region, 
                       data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.niche_exp <- phylostep(lm.niche_exp$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.niche_exp$residuals^2)/sum(null.niche_exp$residuals^2)
1-sum(step.lm.niche_exp$residuals^2)/sum(null.niche_exp$residuals^2)


save(null.rangeD, lm.rangeD, step.lm.rangeD,
     null.range_stab, lm.range_stab, step.lm.range_stab,
     null.range_unf, lm.range_unf, step.lm.range_unf,
     null.range_exp, lm.range_exp, step.lm.range_exp,
     null.nicheD, lm.nicheD, step.lm.nicheD,
     null.niche_stab, lm.niche_stab, step.lm.niche_stab,
     null.niche_unf, lm.niche_unf, step.lm.niche_unf,
     null.niche_exp, lm.niche_exp, step.lm.niche_exp, 
     file = file.path(data_dir,"phylo_trait_models.RData"))


#------------ variable importance ---------

# code adapted from Zurell et al. (2018) JBI


# Helper functions:
R2 <- function(mod0, mod1){
  p <- length(mod1$coefficients)
  n <- length(mod1$fitted)
 
  R2 <- round(1-sum(mod1$residuals^2)/sum(mod0$residuals^2), 3)
  R2adj <- 1 - ((n - 1)/(n - p)) * (1 - R2)
  R2adj <- ifelse(R2adj<0, 0, R2adj)
  return(c(R2=R2, R2adj=R2adj)) 
}

R2glm <- function(model){
  p <- length(model$coefficients)
  n <- length(model$fitted)
  
  R2 <- round(1 - ( model$deviance / model$null.deviance ), 2)
  R2adj <- 1 - ((n - 1)/(n - p)) * (1 - R2)
  R2adj <- ifelse(R2adj<0, 0, R2adj)
  return(c(R2=R2, R2adj=R2adj)) 
}

ETP <- function(x) {eval(parse(text=x))}


# Comparative dataset creation for use in functions - only needed when using package caper (pgls function)
# all_data <- comparative.data(phy=phylo_all, data=traits_all, names.col = "phylo_species", vcv=T, vcv.dim=3)

# number of replicates
nrep=99

covariates <- c("Mass", "Hand.Wing.Index", "Trophic.Level", "Habitat.Density", "Range.Size", "Migration", "Centroid.Latitude", "niche_breadth_zcor", "region")

# names of models
MOD <- c("step.lm.rangeD", "step.lm.range_stab", "step.lm.range_unf", "step.lm.range_exp", "step.lm.nicheD", "step.lm.niche_stab", "step.lm.niche_unf", "step.lm.niche_exp")

output <- list()

for(i in MOD){
  modl <- ETP(i)
  respvar <- strsplit(formula(modl), " ")[[1]][1]
  explvar <- names(modl$coefficients)[2:length(modl$coefficients)]
  # careful - we have a factor level!
  explvar[explvar %in% "regionUS"] <- "region"
  
  # NOTE: if you include quadratic terms, then some more fiddling might be needed here.
  
  # Calculate the R2 with a null model for which we fix the lambda value
  if(modl$optpar>0.000001) {
   # pgls0 <- ETP(paste("pgls(", respvar, "~ 1, data=all_data, lambda=modl$optpar)"))
   pgls0 <- ETP(paste("phylolm(",respvar, "~ 1, data=traits_all, phy=phylo_all, model='lambda')"))
   # pgls1 <- ETP(paste("pgls(", formula(modl), ", data=all_data, lambda=modl$optpar)"))
   pgls1 <- ETP(paste("phylolm(",formula(modl), ", data=traits_all, phy=phylo_all, model='lambda')"))
   r2 <- R2(pgls0, pgls1)
    # Calculate the variable importance with variable randomization and fixed lambda
    rR2 <- rR2a <- data.frame(matrix(NA, nrep, length(explvar), T, list(1:nrep, explvar)))
    for(j in explvar){
      # tempdat <- pgls1$data
      tempdat <- traits_all
      for(k in 1:nrep){	
        # tempdat$data[,j] <- sample(tempdat$data[,j])
        tempdat[,j] <- sample(tempdat[,j])
        # mt <- ETP(paste("pgls(", formula(modl), ", data=tempdat, lambda=modl$optpar)"))
        mt <- ETP(paste("phylolm(",formula(modl), ", data=tempdat, phy=phylo_all, model='lambda')"))
        rR2[k,j] <- R2(pgls0, mt)[1] ; rR2a[k,j] <- R2(pgls0, mt)[2]
        cat(k)
      }
      print(j)
    }
  } else {
    pgls0 <- ETP(paste("glm(", respvar, "~ 1, data=traits_all)"))
    pgls1 <- ETP(paste("glm(", formula(modl), ", data=traits_all)"))
    r2 <- R2glm(pgls1)
    # Calculate the variable importance with variable randomization and fixed lambda
    rR2 <- rR2a <- data.frame(matrix(NA, nrep, length(explvar), T, list(1:nrep, explvar)))
    for(j in explvar){
      tempdat <- pgls1$data
      for(k in 1:nrep){	
        tempdat[,j] <- sample(tempdat[,j])
        mt <- ETP(paste("glm(", formula(modl), ", data=tempdat)"))
        rR2[k,j] <- R2glm(mt)[1] ; rR2a[k,j] <- R2glm(mt)[2]
        cat(k)
      }
      print(j)
    }
  }	
  
  output[[i]] <- list(obs=r2, randR2=rR2, randR2adj=rR2a)
  print(i)
}

# Calculate variable importance (based on mean and SD of simulated R2s)
step.varImp.allmetrics = output
# get variable importances
for (i in 1:length(step.varImp.allmetrics))	{
  tt <- step.varImp.allmetrics[[i]]
  r2 <- tt[[1]][1] - tt$randR2 
  r2a <- tt[[1]][2] - tt$randR2adj
  # Rescale the values (?)
  if (ncol(r2)>1) {
    r2 <- t(apply(r2, 1, function(x) {x=ifelse(x<0,10^-20,x); x/sum(x)} ))
    r2a <- t(apply(r2a, 1, function(x) {x=ifelse(x<0,10^-20,x); x/sum(x)} ))
    # Get the means and sd
    Mr2 <- apply(r2, 2, mean) ; Mr2a <- apply(r2a, 2, mean)
    Sr2 <- apply(r2, 2, sd) ; Sr2a <- apply(r2a, 2, sd) } else {
    r2 = r2/sum(r2)
    r2a=r2a/sum(r2a)
    Mr2 = mean(r2[,1])
    Mr2a = mean(r2a[,1])
    Sr2 = sd(r2[,1])
    Sr2a = sd(r2a[,1])
  }
  step.varImp.allmetrics[[i]]['Mr2'] <- list(Mr2)
  step.varImp.allmetrics[[i]]['Sr2'] <- list(Sr2) 
  step.varImp.allmetrics[[i]]['Mr2a'] <- list(Mr2a)
  step.varImp.allmetrics[[i]]['Sr2a'] <- list(Sr2a)
}

# the metric Mr2 is of main interest (mean R2)
save(step.varImp.allmetrics, file = file.path(data_dir,"VarImp_phylo_trait_models.RData"))

for (i in MOD) {
  print(i)
  print(ETP(paste0("step.varImp.allmetrics","$",i,"$obs[1]")))
  print(ETP(paste0("step.varImp.allmetrics","$",i,"$Mr2")))
}



#-------------- STORE RESULTS ------------------

# Make data.frame to summarise results of trait models
results_TraitAnal_df <- data.frame(matrix(nrow=length(covariates)+3,ncol=length(MOD)*4+1))
names(results_TraitAnal_df) <- c("Trait",paste0(rep(sub("step.lm.","", MOD),each=4),c("_coef","_stderr","_p","_varimp")))
results_TraitAnal_df[,1] <-  c("Intercept",covariates, "R2", "lambda")

# Make vector of null model names
MOD_null <- sub("step.lm","null",MOD)

# loop through models to collect coefficients and variable importance
for (i in 1:length(MOD)) {
  modl <- ETP(MOD[i])
  
  explvar <- names(modl$coefficients)[2:length(modl$coefficients)]
  # careful - we have a factor level!
  explvar[explvar %in% "regionUS"] <- "region"
  
  # store model coefficients
  results_TraitAnal_df[c(1,which(covariates %in% explvar)+1),1+i*4-3] <- round(as.numeric(modl$coefficients),3)
  
  # store R2 of model
  results_TraitAnal_df[nrow(results_TraitAnal_df)-1,1+i*4-3] <- R2(ETP(MOD_null[i]),modl)[[1]]
  
  # store lambda of model
  results_TraitAnal_df[nrow(results_TraitAnal_df),1+i*4-3] <- round(modl$optpar,3)
  
  # store std errors
  results_TraitAnal_df[c(1,which(covariates %in% explvar)+1),1+i*4-2] <- summary(modl)$coefficients[,2]
  
  # store p-values
  results_TraitAnal_df[c(1,which(covariates %in% explvar)+1),1+i*4-1] <- round(as.numeric(summary(modl)$coefficients[,4]),3)
  
  # store variable importance
  varimp <- round(step.varImp.allmetrics[[MOD[i]]]$Mr2,3)
  results_TraitAnal_df[which(covariates %in% explvar)+1,1+i*4] <- ifelse(varimp<0,0,varimp)
}

write.csv(results_TraitAnal_df, file=file.path(data_dir,"results_TraitAnal_df.csv"), row.names=F)


#---------------------------
# trait analyses for separate regions:

traits_all_all <- traits_all
phylo_all_all <- phylo_all

# subset to single region?
reg <- 'US' 
# reg <- 'Europe'
traits_all <- traits_all_all %>%
  filter(region==reg) # 'Europe'
models_filename <- paste0("phylo_trait_models_",reg,".RData")
traitresults_csvname <- paste0("results_TraitAnal_df_",reg,".csv")
traitsimp_filename <- paste0("VarImp_phylo_trait_models_",reg,".RData")

# also remove from phylogeny
phylo_all= drop.tip(phylo_all_all,phylo_all_all$tip.label[!phylo_all_all$tip.label %in% sub(' ','_',traits_all$phylo_species)])


# add rownames - needed for matching phylogenetic information
rownames(traits_all) = sub(' ','_',traits_all$phylo_species)

# Standardise traits
traits_all$Mass <- as.numeric(scale(traits_all$Mass))
traits_all$Hand.Wing.Index <- as.numeric(scale(traits_all$Hand.Wing.Index))
traits_all$Range.Size <- as.numeric(scale(traits_all$Range.Size))
traits_all$Centroid.Latitude <- traits_all$Centroid.Latitude/90
traits_all$Centroid.Longitude <- traits_all$Centroid.Longitude/180
traits_all$niche_breadth_zcor <- as.numeric(scale(traits_all$niche_breadth_zcor))
traits_all$range_shift <- as.numeric(scale(traits_all$range_shift))
traits_all$region <- as.factor(traits_all$region)

# Logit-transform response variables
traits_all$range_D_logit <- logit(traits_all$range_D)
traits_all$range_stability_std_logit <- logit(traits_all$range_stability_std)
traits_all$range_unfilling_std_logit <- logit(traits_all$range_unfilling_std)
traits_all$range_expansion_std_logit <- logit(traits_all$range_expansion_std)
traits_all$niche_D_logit <- logit(traits_all$niche_D)
traits_all$niche_stability_std_logit <- logit(traits_all$niche_stability_std)
traits_all$niche_unfilling_std_logit <- logit(traits_all$niche_unfilling_std)
traits_all$niche_expansion_std_logit <- logit(traits_all$niche_expansion_std)


#-------------- MODELS ---------------

# range_D 
# null model
null.rangeD = phylolm(range_D_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.rangeD = phylolm(range_D_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor, 
                    data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.rangeD <- phylostep(lm.rangeD$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.rangeD$residuals^2)/sum(null.rangeD$residuals^2)
1-sum(step.lm.rangeD$residuals^2)/sum(null.rangeD$residuals^2)
summary(step.lm.rangeD)


# range_stability_std_logit 
# null model
null.range_stab = phylolm(range_stability_std_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.range_stab = phylolm(range_stability_std_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor, 
                        data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.range_stab <- phylostep(lm.range_stab$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.range_stab$residuals^2)/sum(null.range_stab$residuals^2)
1-sum(step.lm.range_stab$residuals^2)/sum(null.range_stab$residuals^2)
summary(step.lm.range_stab)


# range_unfilling_std_logit 
# null model
null.range_unf = phylolm(range_unfilling_std_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.range_unf = phylolm(range_unfilling_std_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor, 
                       data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.range_unf <- phylostep(lm.range_unf$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.range_unf$residuals^2)/sum(null.range_unf$residuals^2)
1-sum(step.lm.range_unf$residuals^2)/sum(null.range_unf$residuals^2)
summary(step.lm.range_unf)


# range_expansion_std_logit 
# null model
null.range_exp = phylolm(range_expansion_std_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.range_exp = phylolm(range_expansion_std_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor, 
                       data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.range_exp <- phylostep(lm.range_exp$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.range_exp$residuals^2)/sum(null.range_exp$residuals^2)
1-sum(step.lm.range_exp$residuals^2)/sum(null.range_exp$residuals^2)
summary(step.lm.range_exp)


# niche_D 
# null model
null.nicheD = phylolm(niche_D_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.nicheD = phylolm(niche_D_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor, 
                    data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.nicheD <- phylostep(lm.nicheD$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.nicheD$residuals^2)/sum(null.nicheD$residuals^2)
1-sum(step.lm.nicheD$residuals^2)/sum(null.nicheD$residuals^2)
summary(step.lm.nicheD)


# niche_stability_std_logit 
# null model
null.niche_stab = phylolm(niche_stability_std_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.niche_stab = phylolm(niche_stability_std_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor, 
                        data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.niche_stab <- phylostep(lm.niche_stab$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.niche_stab$residuals^2)/sum(null.niche_stab$residuals^2)
1-sum(step.lm.niche_stab$residuals^2)/sum(null.niche_stab$residuals^2)
summary(step.lm.niche_stab)


# niche_unfilling_std_logit 
# null model
null.niche_unf = phylolm(niche_unfilling_std_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.niche_unf = phylolm(niche_unfilling_std_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor, 
                       data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.niche_unf <- phylostep(lm.niche_unf$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.niche_unf$residuals^2)/sum(null.niche_unf$residuals^2)
1-sum(step.lm.niche_unf$residuals^2)/sum(null.niche_unf$residuals^2)
summary(step.lm.niche_unf)


# niche_expansion_std_logit 
# null model
null.niche_exp = phylolm(niche_expansion_std_logit ~ 1, data=traits_all, phy=phylo_all, model='lambda')
# full models
lm.niche_exp = phylolm(niche_expansion_std_logit ~ Mass + Hand.Wing.Index + Trophic.Level + Habitat.Density + Range.Size + Migration + Centroid.Latitude + niche_breadth_zcor, 
                       data=traits_all, phy=phylo_all, model='lambda')
# step model
step.lm.niche_exp <- phylostep(lm.niche_exp$formula, data=traits_all, phy=phylo_all, model='lambda')
# Explained deviance
1-sum(lm.niche_exp$residuals^2)/sum(null.niche_exp$residuals^2)
1-sum(step.lm.niche_exp$residuals^2)/sum(null.niche_exp$residuals^2)
summary(step.lm.niche_exp)


save(null.rangeD, lm.rangeD, step.lm.rangeD,
     null.range_stab, lm.range_stab, step.lm.range_stab,
     null.range_unf, lm.range_unf, step.lm.range_unf,
     null.range_exp, lm.range_exp, step.lm.range_exp,
     null.nicheD, lm.nicheD, step.lm.nicheD,
     null.niche_stab, lm.niche_stab, step.lm.niche_stab,
     null.niche_unf, lm.niche_unf, step.lm.niche_unf,
     null.niche_exp, lm.niche_exp, step.lm.niche_exp, 
     file = file.path(data_dir,models_filename))


########
# variable importance 

# number of replicates
nrep=99

covariates <- c("Mass", "Hand.Wing.Index", "Trophic.Level", "Habitat.Density", "Range.Size", "Migration", "Centroid.Latitude", "niche_breadth_zcor")

# names of models
MOD <- c("step.lm.rangeD", "step.lm.range_stab", "step.lm.range_unf", "step.lm.range_exp", "step.lm.nicheD", "step.lm.niche_stab", "step.lm.niche_unf", "step.lm.niche_exp")

output <- list()

for(i in MOD){
  modl <- ETP(i)
  respvar <- strsplit(formula(modl), " ")[[1]][1]
  explvar <- names(modl$coefficients)[2:length(modl$coefficients)]
  # careful - we have a factor level!
  explvar[explvar %in% "regionUS"] <- "region"
  
  # NOTE: if you include quadratic terms, then some more fiddling might be needed here.
  
  # Calculate the R2 with a null model for which we fix the lambda value
  if(modl$optpar>0.000001) {
    # pgls0 <- ETP(paste("pgls(", respvar, "~ 1, data=all_data, lambda=modl$optpar)"))
    pgls0 <- ETP(paste("phylolm(",respvar, "~ 1, data=traits_all, phy=phylo_all, model='lambda')"))
    # pgls1 <- ETP(paste("pgls(", formula(modl), ", data=all_data, lambda=modl$optpar)"))
    pgls1 <- ETP(paste("phylolm(",formula(modl), ", data=traits_all, phy=phylo_all, model='lambda')"))
    r2 <- R2(pgls0, pgls1)
    # Calculate the variable importance with variable randomization and fixed lambda
    rR2 <- rR2a <- data.frame(matrix(NA, nrep, length(explvar), T, list(1:nrep, explvar)))
    for(j in explvar){
      # tempdat <- pgls1$data
      tempdat <- traits_all
      for(k in 1:nrep){	
        # tempdat$data[,j] <- sample(tempdat$data[,j])
        tempdat[,j] <- sample(tempdat[,j])
        # mt <- ETP(paste("pgls(", formula(modl), ", data=tempdat, lambda=modl$optpar)"))
        mt <- ETP(paste("phylolm(",formula(modl), ", data=tempdat, phy=phylo_all, model='lambda')"))
        rR2[k,j] <- R2(pgls0, mt)[1] ; rR2a[k,j] <- R2(pgls0, mt)[2]
        cat(k)
      }
      print(j)
    }
  } else {
    pgls0 <- ETP(paste("glm(", respvar, "~ 1, data=traits_all)"))
    pgls1 <- ETP(paste("glm(", formula(modl), ", data=traits_all)"))
    r2 <- R2glm(pgls1)
    # Calculate the variable importance with variable randomization and fixed lambda
    rR2 <- rR2a <- data.frame(matrix(NA, nrep, length(explvar), T, list(1:nrep, explvar)))
    for(j in explvar){
      tempdat <- pgls1$data
      for(k in 1:nrep){	
        tempdat[,j] <- sample(tempdat[,j])
        mt <- ETP(paste("glm(", formula(modl), ", data=tempdat)"))
        rR2[k,j] <- R2glm(mt)[1] ; rR2a[k,j] <- R2glm(mt)[2]
        cat(k)
      }
      print(j)
    }
  }	
  
  output[[i]] <- list(obs=r2, randR2=rR2, randR2adj=rR2a)
  print(i)
}

# Calculate variable importance (based on mean and SD of simulated R2s)
step.varImp.allmetrics = output
# get variable importances
for (i in 1:length(step.varImp.allmetrics))	{
  tt <- step.varImp.allmetrics[[i]]
  r2 <- tt[[1]][1] - tt$randR2 
  r2a <- tt[[1]][2] - tt$randR2adj
  # Rescale the values (?)
  if (ncol(r2)>1) {
    r2 <- t(apply(r2, 1, function(x) {x=ifelse(x<0,10^-20,x); x/sum(x)} ))
    r2a <- t(apply(r2a, 1, function(x) {x=ifelse(x<0,10^-20,x); x/sum(x)} ))
    # Get the means and sd
    Mr2 <- apply(r2, 2, mean) ; Mr2a <- apply(r2a, 2, mean)
    Sr2 <- apply(r2, 2, sd) ; Sr2a <- apply(r2a, 2, sd) } else {
      r2 = r2/sum(r2)
      r2a=r2a/sum(r2a)
      Mr2 = mean(r2[,1])
      Mr2a = mean(r2a[,1])
      Sr2 = sd(r2[,1])
      Sr2a = sd(r2a[,1])
    }
  step.varImp.allmetrics[[i]]['Mr2'] <- list(Mr2)
  step.varImp.allmetrics[[i]]['Sr2'] <- list(Sr2) 
  step.varImp.allmetrics[[i]]['Mr2a'] <- list(Mr2a)
  step.varImp.allmetrics[[i]]['Sr2a'] <- list(Sr2a)
}

# the metric Mr2 is of main interest (mean R2)
save(step.varImp.allmetrics, file = file.path(data_dir,traitsimp_filename))

for (i in MOD) {
  print(i)
  print(ETP(paste0("step.varImp.allmetrics","$",i,"$obs[1]")))
  print(ETP(paste0("step.varImp.allmetrics","$",i,"$Mr2")))
}



#-------------- STORE RESULTS ------------------

# Make data.frame to summarise results of trait models
results_TraitAnal_df <- data.frame(matrix(nrow=length(covariates)+3,ncol=length(MOD)*4+1))
names(results_TraitAnal_df) <- c("Trait",paste0(rep(sub("step.lm.","", MOD),each=4),c("_coef","_stderr","_p","_varimp")))
results_TraitAnal_df[,1] <-  c("Intercept",covariates, "R2", "lambda")

# Make vector of null model names
MOD_null <- sub("step.lm","null",MOD)

# loop through models to collect coefficients and variable importance
for (i in 1:length(MOD)) {
  modl <- ETP(MOD[i])
  
  explvar <- names(modl$coefficients)[2:length(modl$coefficients)]
  # careful - we have a factor level!
  explvar[explvar %in% "regionUS"] <- "region"
  
  # store model coefficients
  results_TraitAnal_df[c(1,which(covariates %in% explvar)+1),1+i*4-3] <- round(as.numeric(modl$coefficients),3)
  
  # store R2 of model
  results_TraitAnal_df[nrow(results_TraitAnal_df)-1,1+i*4-3] <- R2(ETP(MOD_null[i]),modl)[[1]]
  
  # store lambda of model
  results_TraitAnal_df[nrow(results_TraitAnal_df),1+i*4-3] <- round(modl$optpar,3)
  
  # store std errors
  results_TraitAnal_df[c(1,which(covariates %in% explvar)+1),1+i*4-2] <- summary(modl)$coefficients[,2]
  
  # store p-values
  results_TraitAnal_df[c(1,which(covariates %in% explvar)+1),1+i*4-1] <- round(as.numeric(summary(modl)$coefficients[,4]),3)
  
  # store variable importance
  varimp <- round(step.varImp.allmetrics[[MOD[i]]]$Mr2,3)
  results_TraitAnal_df[which(covariates %in% explvar)+1,1+i*4] <- ifelse(varimp<0,0,varimp)
}

write.csv(results_TraitAnal_df, file=file.path(data_dir,traitresults_csvname), row.names=F)
