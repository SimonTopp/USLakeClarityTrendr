## This file sets up the parameters that are consistent between scripts and then is sourced at the beginning of each markdown


#Parameters
iteration = 'FFS_PolyCorr_noLn_Secchi'

#Set study regions
region <- st_read('in/NLA/NLA_Ecoregions/EcoRegsMerged.shp')

#Set the lake sample subset we're working with
lakeSamp = 'NLA2012_cntr' #'Over10' #"EcoReg2000"

## Source the necessray modelling/analysis functions
source('0_ModellingFunctions.R')

#Non-log transformed model (holdover from data exploration phase)
log = F

##Identify work stage
stage = 5
if(stage > 1){
  #Landsat Correction Models
  load('models/landsat_poly_corrs.Rdata')
}

if(stage > 3){
#Identify the downloaded data and their metadata/filepaths
## Read in lakes sent up to EE
  if(lakeSamp == 'NLA' | lakeSamp == 'NLA2012_cntr'){
    lakes.up <- read_feather('out/lakesNLA2012.feather')
  }else if(lakeSamp == 'EcoReg2000'){  
    lakes.up <- read_feather('out/lakesEcoReg2000.feather')
  }else if(lakeSamp == 'Over10'){
    lakes.up = read_feather('out/Over10sqkm_cntr.feather')
  }
  
  lakesDown <- list.files(paste0('lake_data/',lakeSamp), full.names = T)
  empties <- lakesDown[file.info(lakesDown)[["size"]]==1]
  lakesDown <- lakesDown[!lakesDown %in% empties]
  
  
  if(lakeSamp == 'NLA' | lakeSamp == 'NLA2012_cntr'){
    lake.join <- read_feather('out/NLA2012LakesFull.feather')
  }else if(lakeSamp == 'EcoReg2000'){
    lake.join <- read_feather('out/EcoReg2000LakesFull.feather')
  }else if(lakeSamp == 'Over10'){
    lake.join <- read_feather('out/Over10LakesFull.feather')
  }
  
  ids <- list.files(paste0('lake_data/',lakeSamp)) %>% strsplit(split = '.csv', fixed = T) %>% unlist()
}

if(stage > 4){
  load(paste0('models/',iteration, '.Rdata'))
  Preds.out <- read_feather(paste0('out/TS_Preds/', lakeSamp, '_',iteration,'.feather'))
}