## author: Anita Carolyne Orera
## adapted from: Damien Georges biomod 2 multi species modelling
## date: 2024-01-03


## setup environment ----
## Setup environment ----
setwd('C:/Users/USER/Desktop/geoMountains/multiSpecies/workdir_2')


## load the required packages
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

# ----------------------------------------------------------------- Call the list of current biophysical variables
biophys_var_Cur <-  list.files("../data/biophysical_var/current", pattern = ".asc", full.names = TRUE)


#********************************************************************************************************************************************************************
#----------------------------- RUNNING A SELECTED MODEL
#********************************************************************************************************************************************************************
# 1) -------------------------------------------------------------------- Current biophysical variables
# --- Picking the selected variables as rasters (Refer to single species models)
pop <- raster(biophys_var_Cur[[5]]) #This raster file had a different X estent. All other rasters would have to be modified to fit its extent
decNDVI <- raster(biophys_var_Cur[[1]])
julNDVI <- raster(biophys_var_Cur[[2]])

# ******  Reset the extent 
# Step 1 - Crop the rasters to fit in the X and Y limits of the population dataset
decNDVI <- crop(decNDVI, extent(pop))
julNDVI <- crop(julNDVI, extent(pop))
# Step 2 - Resample the rasters to fit in the X and Y dimensions of the population dataset
decNdvi <- resample(decNDVI, pop, method = "bilinear")
julNdvi <- resample(julNDVI, pop, method = "bilinear")
# Step 3 - Set the extent of the current stack to match that of the projected stack
extent(decNdvi) <- extent(pop)
extent(julNdvi) <- extent(pop)

#--- Stack the selected variables
biophys_KE_current <- stack(
                                c(
                                    ndvi_2 = decNdvi,
                                    ndvi_1 = julNdvi,
                                    pop_pul = pop
                                ),
                                RAT = FALSE
                           )


# 2) ------------------------------------------------------------------- 2050 biophysical variables (rcp 2.6)
## Load 2050 biophysical variables (RCP 2.6)
biophys_KE_2050_BC26 <- stack(
                              c(
                                ndvi_2 = "../data/biophysical_var/2050_26/meanndviextent_conservancies_final.asc",
                                pop_pul = "../data/biophysical_var/2050_26/ken_pop_pul_2012_conservancies_final.asc",                          
                                ndvi_1 = "../data/biophysical_var/2050_26/stdndvi_extent_conservancies_final.asc"
                               ),
                               RAT = FALSE
                             )


# 3) ------------------------------------------------------------------- 2050 biophysical variables (rcp 8.5)
## Load 2050 biophysical variables (RCP 8.5)
biophys_KE_2050_BC85 <- stack(
                              c(
                                ndvi_2 = "../data/biophysical_var/2050_85/meanndviextent_conservancies_final.asc",
                                pop_pul = "../data/biophysical_var/2050_85/ken_pop_pul_2012_conservancies_final.asc",                          
                                ndvi_1 = "../data/biophysical_var/2050_85/stdndvi_extent_conservancies_final.asc"
                               )
                             )


# 4) ------------------------------------------------------------------- Build a species modelling wrapper ----
biomod2_wrapper <- function(sp){
                                    cat("\n> species : ", sp)
                                    
                                    ## Get occurrences points
                                    sp_dat <- data[data$species == sp, ]
                                    ## Format the data
                                    sp_format <- BIOMOD_FormatingData(
                                                                        resp.var = rep(1, nrow(sp_dat)), 
                                                                        expl.var = biophys_KE_current,
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
                                                                                metric.select.thresh = 0.75,
                                                                                metric.eval = c('TSS','ROC'),
                                                                                var.import = 0,
                                                                                em.algo = c('EMmean', 'EMcv', 'EMca') #You can add 'EMwmean' but in our case it produced null results
                                                                           )
                                    
                                    ## Do projections
                                    proj_scen <- c("current", "2050_BC26", "2050_BC85")
                                    for(scen in proj_scen){
                                                            cat("\n> projections of ", scen)
                                                            ## single model projections
                                                            sp_proj <- BIOMOD_Projection(
                                                                                            bm.mod = sp_model,
                                                                                            new.env = get(paste0("biophys_KE_", scen)),
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
#---- The step for alpha diversity map production was subjectively completed with mean.
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

##---------------------------- 2050 (RCP 2.6) conditons
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
physMeanPlot <- levelplot(
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
physMeanPlot 