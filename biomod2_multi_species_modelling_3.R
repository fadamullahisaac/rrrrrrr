## author: Anita Carolyne
## adapted from: Damien Georges biomod 2 multi species modelling
## date: 2024-01-05


## Setup environment ----
setwd('C:/Users/USER/Desktop/geoMountains/multiSpecies/workdir_3')


## Load the required packages
library(biomod2)
library(raster)
library(rasterVis)
library(gridExtra)
library(reshape2)


#********************************************************************************************************************************************************************
#-----------------------------READING THE RELEVANT DATA
#********************************************************************************************************************************************************************
## Read data ----
## ----------------------------------------------------------------- Species occurences data
data <- read.csv('../data/species.csv', stringsAsFactors = FALSE)
head(data) #--- Displays the first six entries
table(data$species) # Displays the entry count per species - recorded instances
spp_to_model <- unique(data$species) #--- Displays the unique species names - names of all the species we want to model

# ----------------------------------------------------------------- Define an output path for R outputs
outPath <- "../outputs"
# ----------------------------------------------------------------- Call the list of; 
# ******** Current biophysical variables
biophys_var_Cur <-  list.files("../data/biophysical_var/current", pattern = ".asc", full.names = TRUE)
# ------A. Store the biophysical rasters as variables ------ 
pop <- raster(biophys_var_Cur[[5]]) #This raster file had a different X estent. Other rasters would have to be modified to fit its extent
meanNDVI <- raster(biophys_var_Cur[[6]])
# ******B.  Reset the extent 
meanNDVI <- crop(meanNDVI, extent(pop))
# ******C.  Perform resampling
meanNdvi <- resample(meanNDVI, pop, method = "bilinear")
# ******D.  Perform matches
extent(meanNdvi) <- extent(pop)


# 1) CURRENT SCENARIO
# a) -------------------------------------------------------------------- Current climatic variables
bioclim <- stack(
                                c(
                                    bio4 = '../data/current_1960_1990_ASC/worldclim_KE_bio4.asc',
                                    bio3 = '../data/current_1960_1990_ASC/worldclim_KE_bio3.asc',
                                    bio2  = '../data/current_1960_1990_ASC/worldclim_KE_bio2.asc',
                                    bio6  = '../data/current_1960_1990_ASC/worldclim_KE_bio6.asc', 
                                    bio14 = '../data/current_1960_1990_ASC/worldclim_KE_bio14.asc',
                                    bio7 = '../data/current_1960_1990_ASC/worldclim_KE_bio7.asc',
                                    bio18 = '../data/current_1960_1990_ASC/worldclim_KE_bio18.asc',
                                    bio15 = '../data/current_1960_1990_ASC/worldclim_KE_bio15.asc' 
                                ),
                                RAT = FALSE
                )
# b) -------------------------------------------------------------------- Current biophysical variables
# ------ Stack the corrected biophysical ASC files
# -------------- Biophysical variables
bioPhys <-  stack(
                  c(
                        pop_pul = pop,
                        mean_ndvi = meanNdvi
                   ),
                    RAT = FALSE
                 )

# ******)  Call to print
# ----- Note that our bioPhys stack now has a different extent to that of the bioclim stack
bioclim
bioPhys
# ----- To correct this, we set the bioclim stack to have the same extent and resolution as the bioPhys stack

# *****) Check for consistency in dimensions (extent & resolution)
# Step 1 - Perform a dimensions check
dim(bioclim)
dim(bioPhys)
# Step 2 - Perform a resolution check
res(bioclim)
res(bioPhys)
# Step 3 - check the extent by min & max values
bioclim@extent
extent(bioPhys)

# *****) Reset the extent & resolution
# Step 1 - Crop the current bioclimatic raster stack to fit in the X and Y limits of the biophysical dataset
bioclim_cropped <- crop(bioclim, extent(bioPhys))
# Step 2 - Set the resolution of the current cropped stack to match that of the projected stack
res(bioclim_cropped) <- res(bioPhys)
# Step 3 - Set the extent of the current stack to match that of the biophysical stack
extent(bioclim_cropped) <- extent(bioPhys)
# Step 4 - Copy the values from the cropped stack (retention of data)
bioclim_cropped[] <- bioclim_cropped[]
bioclim_matched <- bioclim_cropped

# *****) Confirm for consistency in extent & resolution
# Step 1 - Resolution check
res(bioclim_matched)
res(bioPhys)
# Step 2 - Extent check
bioclim_matched@extent
extent(bioPhys)
# Step 3 - Resample the stack to fit in the X and Y dimensions of the bioPhys stack
bioClim <- resample(bioclim_matched, bioPhys, method = "bilinear")


# *****)  Stack both the bioclimatic and biophysical variables
bioClim_bioPhys <- stack(bioClim, bioPhys)


# c) -------------------------------------------------------------------- Combined data in use 
# --- By first using the names function to allow subsetting in the correct named order
names(bioClim_bioPhys)
bioClim_bioPhys_current <- stack(
                                    c(
                                      bio4 = bioClim_bioPhys[[1]],
                                      bio3  = bioClim_bioPhys[[2]],
                                      bio2 = bioClim_bioPhys[[3]],
                                      pop = bioClim_bioPhys[[9]],
                                      bio6 = bioClim_bioPhys[[4]],
                                      bio14 = bioClim_bioPhys[[5]],                                      
                                      bio7  = bioClim_bioPhys[[6]],                               
                                      bio15 = bioClim_bioPhys[[8]],
                                      bio18 = bioClim_bioPhys[[7]],
                                      ndvi  = bioClim_bioPhys[[10]]
                                    )
                                  )
# ------ Call to print
bioClim_bioPhys_current


# 2) FUTURE SCENARIOS
# 2.1) ------------------------------------------------------------------- 2050 variables (rcp 2.6)
## Load 2050 bioclim variables (RCP 2.6)
bioclim_KE_2050_BC26 <- stack(
                              c(
                                bio4 = '../data/worldclim_KE/worldclim_KE_2050_BC26_Bio4.asc',
                                bio3 = '../data/worldclim_KE/worldclim_KE_2050_BC26_Bio3.asc',
                                bio2  = '../data/worldclim_KE/worldclim_KE_2050_BC26_Bio2.asc', 
                                bio6  = '../data/worldclim_KE/worldclim_KE_2050_BC26_Bio6.asc',
                                bio14 = '../data/worldclim_KE/worldclim_KE_2050_BC26_Bio14.asc',                   
                                bio7 = '../data/worldclim_KE/worldclim_KE_2050_BC26_Bio7.asc',
                                bio18 = '../data/worldclim_KE/worldclim_KE_2050_BC26_Bio18.asc', 
                                bio15 = '../data/worldclim_KE/worldclim_KE_2050_BC26_Bio15.asc'
                               ),
                               RAT = FALSE
                             )
## Load 2050 bioPhys variables (RCP 2.6)
biophys_KE_2050_BC26 <- stack(
                              c(
                                pop = "../data/biophysical_var/2050_26/ken_pop_pul_2012_conservancies_final.asc",
                                ndvi = "../data/biophysical_var/2050_26/meanndviextent_conservancies_final.asc"
                               ),
                               RAT = FALSE
                             )
## Stack the 2050 biophysical and bioclimatic variables (RCP 2.6)
bioClim_bioPhys_2050_BC26 <- stack(bioclim_KE_2050_BC26, biophys_KE_2050_BC26)

# 2.2) ------------------------------------------------------------------- 2050 variables (rcp 8.5)
## Load 2050 bioclim variables (RCP 8.5)
bioclim_KE_2050_BC85 <- stack(
                              c(
                                bio4 = '../data/worldclim_KE/worldclim_KE_2050_BC85_Bio4.asc',
                                bio3 = '../data/worldclim_KE/worldclim_KE_2050_BC85_Bio3.asc',
                                bio2  = '../data/worldclim_KE/worldclim_KE_2050_BC85_Bio2.asc',
                                bio6  = '../data/worldclim_KE/worldclim_KE_2050_BC85_Bio6.asc', 
                                bio14 = '../data/worldclim_KE/worldclim_KE_2050_BC85_Bio14.asc',
                                bio7 = '../data/worldclim_KE/worldclim_KE_2050_BC85_Bio7.asc',
                                bio18 = '../data/worldclim_KE/worldclim_KE_2050_BC85_Bio18.asc',
                                bio15 = '../data/worldclim_KE/worldclim_KE_2050_BC85_Bio15.asc'
                               ),
                               RAT = FALSE
                             )
## Load 2050 biophys variables (RCP 8.5)
biophys_KE_2050_BC85 <- stack(
                              c(
                                pop = "../data/biophysical_var/2050_85/ken_pop_pul_2012_conservancies_final.asc",
                                ndvi = "../data/biophysical_var/2050_85/meanndviextent_conservancies_final.asc" 
                               )
                             )
## Stack the 2050 biophysical and bioclimatic variables (RCP 8.5)
bioClim_bioPhys_2050_BC85 <- stack(bioclim_KE_2050_BC85, biophys_KE_2050_BC85)


# 3) ------------------------------------------------------------------- Build a species modelling wrapper ----
biomod2_wrapper <- function(sp){
                                    cat("\n> species : ", sp)
                                    
                                    ## Get occurrences points
                                    sp_dat <- data[data$species == sp, ]
                                    ## Format the data
                                    sp_format <- BIOMOD_FormatingData(
                                                                        resp.var = rep(1, nrow(sp_dat)), 
                                                                        expl.var = bioClim_bioPhys_current,
                                                                        resp.xy = sp_dat[, c("Long", "Lat")],
                                                                        resp.name = sp,
                                                                        PA.strategy = "random", 
                                                                        PA.nb.rep = 2, 
                                                                        PA.nb.absences = 500
                                                                     )
                                    ## Print the formatting summary
                                    sp_format
                                    
                                    ## Save image of input data summary
                                    if(!exists(sp)) dir.create(sp)
                                    pdf(paste(sp, "/", sp ,"_data_formated.pdf", sep="" ))
                                    try(plot(sp_format))
                                    dev.off()
                                    
                                    ## Define models options
                                    sp_opt <- BIOMOD_ModelingOptions()
                                    
                                    ## Model the species
                                    sp_model <- BIOMOD_Modeling( 
                                                                    bm.format = sp_format, 
                                                                    modeling.id = "demoA",
                                                                    models = c('RF', 'GBM', 'FDA'), 
                                                                    bm.options = sp_opt, 
                                                                    CV.nb.rep = 2, 
                                                                    CV.perc = 0.7, 
                                                                    CV.do.full.models = FALSE,
                                                                    weights = NULL,
                                                                    metric.eval = c('TSS', 'ROC'),
                                                                    var.import = 3,
                                                                    scale.models = FALSE
                                                               )
                                    
                                    ## Save some graphical outputs
                                    #### models scores
                                    pdf(paste0(sp, "/", sp , "_models_scores.pdf"))
                                    try(gg1 <- bm_PlotEvalMean(sp_model, metric.eval = c("TSS", "ROC"), group.by = 'models', do.plot = TRUE))
                                    try(gg2 <- bm_PlotEvalMean(sp_model, metric.eval = c("TSS", "ROC"), group.by = 'data_set', do.plot = TRUE))
                                    try(gg3 <- bm_PlotEvalMean(sp_model, metric.eval = c("TSS", "ROC"), group.by = 'cv_run', do.plot = TRUE))
                                    try(grid.arrange(gg1, gg2, gg3))
                                    dev.off()
                                    
                                    ## Build ensemble models
                                    sp_ens_model <- BIOMOD_EnsembleModeling(
                                                                                bm.mod = sp_model,
                                                                                em.by = 'all',
                                                                                metric.select = 'TSS',
                                                                                metric.select.thresh = 0.75, #A lower threshold value results in missing arguments returned to min max for the selected case
                                                                                metric.eval = c('TSS','ROC'),
                                                                                var.import = 0,
                                                                                em.algo = c('EMmean', 'EMcv', 'EMca', 'EMwmean')
                                                                           )
                                    
                                    ## Do projections
                                    proj_scen <- c("current", "2050_BC26", "2050_BC85")
                                    for(scen in proj_scen){
                                                            cat("\n> projections of ", scen)
                                                            ## Single model projections
                                                            sp_proj <- BIOMOD_Projection(
                                                                                            bm.mod = sp_model,
                                                                                            new.env = get(paste0("bioClim_bioPhys_", scen)),
                                                                                            proj.name = scen,
                                                                                            models.chosen = 'all',
                                                                                            metric.binary = "TSS",
                                                                                            metric.filter = NULL,
                                                                                            compress = TRUE,
                                                                                            build.clamping.mask = FALSE,
                                                                                            do.stack = FALSE,
                                                                                            output.format = ".img"
                                                                                        )
                                                            ## Ensemble model projections
                                                            sp_ens_proj <- BIOMOD_EnsembleForecasting(
                                                                                                        bm.em = sp_ens_model,
                                                                                                        bm.proj = sp_proj,
                                                                                                        metric.binary = "TSS",
                                                                                                        compress = TRUE,
                                                                                                        do.stack = FALSE,
                                                                                                        output.format = ".img"
                                                                                                     )
                                                        }
                                    return(paste0(sp," modelling completed !"))
                                }


# 5) ------------------------------------------------------------------- # Launch the species modelling wrapper over species list ----
## sequential computation
for (sp in spp_to_model){
                            biomod2_wrapper(sp)
                        }
## or with a lapply function in sequential model
## all_species_bm <- lapply(spp_to_model, biomod2_wrapper)
                            

# 6) ------------------------------------------------------------------- # Produce alpha-diversity maps ------------------------------------------
##---------------------------- Current conditions
### -------------------------------- Load binary projections
#---- Using weighted mean returned NA values in all the scenarios, whereas CV lacks a TSS binary image.
#---- The step for alpha diversity map production was thus completed with mean and ca.
#-------------------------------------------------------------------------------------------------------------------------------------------------
# 6.1) ------------------------------------------------------------------- We define the path to the MEAN in binary format
#-------------------------------------------------------------------------------------------------------------------------------------------------
## ---------------------------- Current conditions
f_em_mean_bin_current <- paste0(
                                    spp_to_model,
                                    "/proj_current/individual_projections/", 
                                    spp_to_model, "_EMmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img"
                                )
### sum all projections
if(length(f_em_mean_bin_current) >= 2){
                                            ## initialisation
                                            taxo_alpha_div_current <- raster(f_em_mean_bin_current[1]) 
                                            for(f in f_em_mean_bin_current[-1]){
                                                                                    taxo_alpha_div_current <- taxo_alpha_div_current + raster(f)
                                                                                }                           
                                      }

##---------------------------- 2050 (RCP 2.6) conditions
### load binaries projections
f_em_mean_bin_2050_26 <- paste0(
                                    spp_to_model,
                                    "/proj_2050_BC26/individual_projections/", 
                                    spp_to_model, "_EMmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img"
                                )
### sum all projections
if(length(f_em_mean_bin_2050_26) >= 2){
                                        ## initialisation
                                        taxo_alpha_div_2050_26 <- raster(f_em_mean_bin_2050_26[1]) 
                                        for(f in f_em_mean_bin_2050_26[-1]){
                                                                                taxo_alpha_div_2050_26 <- taxo_alpha_div_2050_26 + raster(f)
                                                                            }
                                      }

##---------------------------- 2050 (RCP 8.5) conditions
### load binaries projections
f_em_mean_bin_2050_85 <- paste0(
                                    spp_to_model,
                                    "/proj_2050_BC85/individual_projections/", 
                                    spp_to_model, "_EMmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img"
                                )
### sum all projections
if(length(f_em_mean_bin_2050_85) >= 2){
                                        ## initialisation
                                        taxo_alpha_div_2050_85 <- raster(f_em_mean_bin_2050_85[1]) 
                                        for(f in f_em_mean_bin_2050_85){
                                                                            taxo_alpha_div_2050_85 <- taxo_alpha_div_2050_85 + raster(f)
                                                                        }               
                                       }


# 7.1) ------------------------------------------------------------------- Perform a full check on the resolution & dimensions
taxo_alpha_div_current
taxo_alpha_div_2050_26
taxo_alpha_div_2050_85
# Resample the current raster to fit in the X and Y dimensions of the 8.5 future dataset (had a higher maximum for values)
species_current <- resample(taxo_alpha_div_current, taxo_alpha_div_2050_85, method = "bilinear")
# Perform a check
species_current


# 8.1) ------------------------------------------------------------------- Plot the alpha-div maps
comboMeanPlot <- levelplot(
                                stack(
                                        c(
                                            current_scenario = species_current, 
                                            in_2050_RCP_2.6 = taxo_alpha_div_2050_26, 
                                            in_2050_RCP_8.5 = taxo_alpha_div_2050_85
                                         )
                                    ),
                                    main = expression(paste("Acacia & Opuntia ", alpha, "-diversity")),
                                    par.settings = BuRdTheme
                            )
comboMeanPlot 


#-------------------------------------------------------------------------------------------------------------------------------------------------
# 6.2) We define the path to the COMNITTEE AVERAGING in binary format
#-------------------------------------------------------------------------------------------------------------------------------------------------
# A "_ca" suffix was added to differentiate from the mean product
f_em_ca_bin_current <- paste0(
                                spp_to_model,
                                "/proj_current/individual_projections/", 
                                spp_to_model, "_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img"
                             )
### sum all projections
if(length(f_em_ca_bin_current) >= 2){
                                            ## initialisation
                                            taxo_alpha_div_current_ca <- raster(f_em_ca_bin_current[1]) 
                                            for(f in f_em_ca_bin_current[-1]){
                                                                                taxo_alpha_div_current_ca <- taxo_alpha_div_current_ca + raster(f)
                                                                             }                           
                                    }

##---------------------------- 2050 (RCP 2.6) conditions
### load binaries projections
f_em_ca_bin_2050_26 <- paste0(
                                spp_to_model,
                                "/proj_2050_BC26/individual_projections/", 
                                spp_to_model, "_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img"
                             )
### sum all projections
if(length(f_em_ca_bin_2050_26) >= 2){
                                        ## initialisation
                                        taxo_alpha_div_2050_26_ca <- raster(f_em_ca_bin_2050_26[1]) 
                                        for(f in f_em_ca_bin_2050_26[-1]){
                                                                            taxo_alpha_div_2050_26_ca <- taxo_alpha_div_2050_26_ca + raster(f)
                                                                         }
                                    }

##---------------------------- 2050 (RCP 8.5) conditons
### load binaries projections
f_em_ca_bin_2050_85 <- paste0(
                                spp_to_model,
                                "/proj_2050_BC85/individual_projections/", 
                                spp_to_model, "_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.img"
                             )
### sum all projections
if(length(f_em_ca_bin_2050_85) >= 2){
                                        ## initialisation
                                        taxo_alpha_div_2050_85_ca <- raster(f_em_ca_bin_2050_85[1]) 
                                        for(f in f_em_ca_bin_2050_85){
                                                                        taxo_alpha_div_2050_85_ca <- taxo_alpha_div_2050_85_ca + raster(f)
                                                                     }               
                                    }


# 7.2) ------------------------------------------------------------------- Perform a full check on the resolution & dimensions
taxo_alpha_div_current_ca
taxo_alpha_div_2050_26_ca
taxo_alpha_div_2050_85_ca
# Resample the current raster to fit in the X and Y dimensions of the 8.5 future dataset (had a higher maximum for values)
species_current_ca <- resample(taxo_alpha_div_current_ca, taxo_alpha_div_2050_85_ca, method = "bilinear")
# Perform a check
species_current_ca


# 8.2) ------------------------------------------------------------------- Plot the alpha-div maps
comboCAplot <- levelplot(
                                stack(
                                        c(
                                            current_scenario = species_current_ca, 
                                            in_2050_RCP_2.6 = taxo_alpha_div_2050_26_ca, 
                                            in_2050_RCP_8.5 = taxo_alpha_div_2050_85_ca
                                         )
                                     ),
                                    main = expression(paste("Acacia & Opuntia ", alpha, "-diversity")),
                                    par.settings = BuRdTheme
                          )
comboCAplot 

                