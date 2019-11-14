## This file sets up the parameters that are consistent between scripts and then is sourced at the beginning of each markdown
library(tidyverse)
library(feather)
library(viridis)
library(knitr)
library(sf)
library(rgdal)
library(maps)
library(magrittr)
library(mlbench)
library(caret)
library(doParallel)
library(furrr)
library(lubridate)
library(groupdata2)
library(CAST)
library(mapview)
library(onehot)
library(Metrics)
library(kableExtra)
library(ggpmisc)
library(reticulate)
library(googledrive)
library(hydrolinks)
library(rkt)
library(xgboost)
library(broom)
library(furrr)
library(plotly)
library(corrgram)
library(boot)
library(trend)
library(scales)
library(gridExtra)

#Parameters
iteration = 'FFS_PolyCorr_noLn_Secchi'

#Set study regions
region <- st_read('in/NLA/NLA_Ecoregions/EcoRegsMerged.shp')

#Set the lake sample subset we're working with
lakeSamp ='NLA2012_cntr' #'NLA2012_cntr' #'Over10' #"EcoReg2000_cntr"

## Source the necessray modelling/analysis functions
source('0_ModellingFunctions.R')

#Non-log transformed model (holdover from data exploration phase)
log = F

##Identify work stage
stage = 6

if(stage > 1){
  
  #Landsat Correction Models
  load('models/landsat_poly_corrs.Rdata')
  
  srMunged <- read_feather('out/srMunged.feather') %>%
    filter(parameter == 'secchi',
           value <= 10,
           !is.na(CatAreaSqKm)) %>% #100 obs that didn't match up with LakeCat
    mutate(COMID = as.character(COMID),
           pctForest2006 = PctDecid2006Cat + PctConif2006Cat + PctMxFst2006Cat,
           pctUrban2006 = PctUrbMd2006Cat + PctUrbHi2006Cat,
           pctWetland2006 = PctWdWet2006Cat + PctHbWet2006Cat,
           areasqkm = round(areasqkm, 1),
           meandused = round(meandused, 1),
           inStreamCat = as.numeric(inStreamCat),
           sat = factor(sat, levels = c(5,7,8), labels = c('l5', 'l7', 'l8'))) 
  
  refCorrect <- srMunged %>% select(UniqueID, blue, green, red, nir, sat) %>%
    gather(blue:nir, key = 'band' , value = 'value') %>%
    spread(sat, value) %>%
    group_by(band) %>%
    nest() %>%
    left_join(funcs.8) %>% #From 1_nhd_join_and_munge
    left_join(funcs.5) %>%
    mutate(pred8 = map2(lm8, data, predict),
           pred5 = map2(lm5, data, predict)) %>%
    select(-c(lm8,lm5)) %>%
    unnest(data, pred8, pred5) %>%
    select(-c(l8,l5)) %>%
    rename(l8 = pred8, l5 = pred5) %>% gather(l5,l7,l8, key = 'sat', value = 'value') %>%
    spread(band, value) %>%
    na.omit()
  
  srMunged <- srMunged %>% select(-c(blue, green, red, nir, sat)) %>% left_join(refCorrect) %>%
    mutate(NR = nir/red,
           BG = blue/green,
           dWL = fui.hue(red, green, blue))
  
  srMunged[srMunged == -9998] = NA ##Standardize missing values
  ##See if removing the ~2k obs with missing depth helps, it doesn't
  #srMunged <- srMunged %>% filter(!is.na(meandused))
  
  if(log == T){srMunged <- srMunged %>% mutate(value = log(value))}
}

if(stage > 2){
  load(paste0('models/',iteration, '.Rdata'))
  output <- read_feather(paste0('out/outputs/',iteration, '.feather'))
}

if(stage > 4){
#Identify the downloaded data and their metadata/filepaths
## Read in lakes sent up to EE
  if(lakeSamp == 'NLA' | lakeSamp == 'NLA2012_cntr'){
    lakes.up <- read_feather('out/lakesNLA2012.feather')
  }else if(lakeSamp == 'EcoReg2000' | lakeSamp == 'EcoReg2000_cntr'){  
    lakes.up <- read_feather('out/lakesEcoReg2000.feather')
  }else if(lakeSamp == 'Over10'){
    lakes.up = read_feather('out/Over10sqkm_cntr.feather')
  }
  
  lakesDown <- list.files(paste0('lake_data/',lakeSamp), full.names = T)
  empties <- lakesDown[file.info(lakesDown)[["size"]]==1]
  lakesDown <- lakesDown[!lakesDown %in% empties]
  
  
  if(lakeSamp == 'NLA' | lakeSamp == 'NLA2012_cntr'){
    lake.join <- read_feather('out/NLA2012LakesFull.feather')
  }else if(lakeSamp == 'EcoReg2000' | lakeSamp == 'EcoReg2000_cntr'){
    lake.join <- read_feather('out/EcoReg2000LakesFull.feather')
  }else if(lakeSamp == 'Over10'){
    lake.join <- read_feather('out/Over10LakesFull.feather')
  }
  
  ids <- list.files(paste0('lake_data/',lakeSamp)) %>% strsplit(split = '.csv', fixed = T) %>% unlist()
}

if(stage > 5){
  Preds.out <- read_feather(paste0('out/TS_Preds/', lakeSamp, '_',iteration,'.feather'))
}
