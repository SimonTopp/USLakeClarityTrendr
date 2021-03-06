
## Function for mapping over evaluation lakes data and generating predictions
EvalPreds <- function(id, paths, lakesUp, log, model, features, lakeSamp){

  # Really buggy but no serious ones, set it up with tryCatch so it'll return an empty data frame if it fails
  predicted <- tibble(system.index = NA, COMID = id, date = NA, year = NA, month = NA, value = NA, lakeSamp = lakeSamp)
  try({
    # df <- read.csv(grep(paste0('/',id,'.csv'),paths, value = T), stringsAsFactors = F) %>%
    # mutate(COMID = as.character(COMID),
    #        sat = factor(sat, levels = c(5,7,8), labels = c('l5','l7', 'l8'))) %>%
    #   filter(!is.na(blue),
    #          dswe == 1,
    #          dswe_sd < .4,
    #          hillshadow == 1,
    #          cScore == 0,
    #          Cloud_Cover < 50,
    #          pixelCount > 5) %>%
    #   filter_at(vars(red, blue, green, nir, swir1, swir2), all_vars(.>0 & .<2000)) %>%
    #   gather(blue, green, red, nir, key = 'band' , value = 'value') %>%
    #   spread(sat, value) %>%
    #   group_by(band) %>%
    #   nest() %>%
    #   left_join(funcs.5) %>%
    #   mutate(pred5 = map2(lm5, data, predict)) %>%
    #   select(-lm5) %>%
    #   unnest(data, pred5) %>%
    #   select(-l5) %>%
    #   rename(l5 = pred5) %>% gather(l5,l7, key = 'sat', value = 'value') %>%
    #   spread(band, value) %>%
    #   filter(!is.na(blue)) %>%
    #   inner_join(lakesUp %>% rename(date.field = date), by = 'COMID') %>%
    #   distinct(.keep_all = T) %>%
    #   mutate(NR = nir/red,
    #          BG = blue/green,
    #          dWL = fui.hue(red, green, blue),
    #          #date = ymd_hms(date),
    #          year = year(date),
    #          month = as.numeric(month(date)),
    #          pctForest2006 = PctDecid2006Cat + PctConif2006Cat + PctMxFst2006Cat,
    #          pctUrban2006 = PctUrbMd2006Cat + PctUrbHi2006Cat,
    #          pctWetland2006 = PctWdWet2006Cat + PctHbWet2006Cat,
    #          areasqkm = round(areasqkm, 1),
    #          meandused = round(meandused, 1)) %>%
    #   filter_at(vars(blue,green,red,nir,swir1,swir2),all_vars(.>0 & .< 2000))


    df <- read.csv(grep(paste0('/',id,'.csv'),paths, value = T), stringsAsFactors = F) %>%
      mutate(COMID = as.character(COMID),
             year = year(date),
             UniqueID = row_number(),
             sat = factor(sat, levels = c(5,7,8), labels = c('l5','l7','l8'))) %>%
      filter(!is.na(blue),
             dswe == 1,
             dswe_sd < .4,
             cScore == 0) %>%
      gather(blue, green, red, nir, key = 'band' , value = 'value') %>%
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
      filter(!is.na(blue)) %>%
      inner_join(lakesUp %>% mutate(COMID = as.character(COMID))) %>%
      distinct(.keep_all = T) %>%
      mutate(NR = nir/red,
             BG = blue/green,
             dWL = fui.hue(red, green, blue),
             #date = ymd_hms(date),
             month = as.numeric(month(date)),
             pctForest2006 = PctDecid2006Cat + PctConif2006Cat + PctMxFst2006Cat,
             pctUrban2006 = PctUrbMd2006Cat + PctUrbHi2006Cat,
             pctWetland2006 = PctWdWet2006Cat + PctHbWet2006Cat,
             areasqkm = round(areasqkm, 1),
             meandused = round(meandused, 1)) %>%
      filter(pixelCount > 5) %>%
      filter_at(vars(blue,green,red,nir,swir1,swir2),all_vars(.>0 & .< 2000))


    if('AOD' %in% features){
      if(lakeSamp == 'NLA'){
      lut <- read_feather('out/aodLUT_NLA2012.feather')
      df <- df %>%
        left_join(lut, by = c('COMID', 'date'))
      }else if(lakeSamp == 'EcoReg2000'){
        lut <- read_feather('out/aodLUT_EcoReg2000.feather')
      df <- df %>%
        left_join(lut)
      }
    }
    
    encoder <- onehot(df %>% select(features))
    lake.input <- predict(encoder, df %>% select(features))

    if(log == T){
      value = tibble(value = exp(predict(model, lake.input)))
    }else{
      value = tibble(value = predict(model, lake.input))
    }
    
    predicted <- df %>%
      select(system.index, COMID, date, year, month, region) %>%
      bind_cols(value) %>%
      mutate(lakeSamp = lakeSamp)

  })
  return(predicted)
}

#Function for making figures over space and time
spaceTimeFigs <- function(con = 'Secchi', abs = T){
  if(abs == T){
    ggplot(errorSum %>% filter(Metric == 'mae'|Metric == 'bias'), aes(x = fct_reorder(quantLabs, order), y = Error, color = Metric, group = interaction(Metric, Variable, log), linetype = log )) + 
      geom_point() + 
      geom_line() + 
      facet_wrap(~Variable, scales = 'free')  +  
      theme_bw() +
      #scale_y_continuous(labels = scales::percent) +
      theme(axis.text.x = element_text(angle = 90, vjust = .5),
            legend.position = 'bottom') +
      labs(x = 'Quantile', y = 'Error (m)', title = paste0(con,' Absolute Errors'))
  }else{
    ggplot(errorSum %>% filter(Metric == 'smape'|Metric == 'p.bias'), aes(x = fct_reorder(quantLabs, order), y = Error, color = Metric, group = interaction(Metric, Variable, log), linetype = log)) + 
      geom_point() + 
      geom_line() + 
      facet_wrap(~Variable, scales = 'free')  +  
      theme_bw() +
      scale_y_continuous(labels = scales::percent) +
      theme(axis.text.x = element_text(angle = 90, vjust = .5),
            legend.position = 'bottom') +
      labs(x = 'Quantile', y = 'Error (m)', title = paste0(con, ' Relative Errors'))
  }
}



## Function for pulling Mera AOD values
pullMerra <- function(path = 'in/MERRA2_data2', date, lat, long){
  files = list.files(path, full.names = T)
  yearMonth <- ifelse(month(date) < 10, paste0(year(date),'0', month(date)),paste0(year(date), month(date)))
  file <- grep(x = files, pattern = yearMonth, value = T)
  brick <- raster::brick(file, varname = 'AODANA')
  AOD <- raster::extract(brick ,SpatialPoints(coords = cbind(long, lat)))[1]
  AOD <- round(AOD, 3)
  return(AOD)
}


fui.hue <- function(R, G, B) {
  
  # Convert R,G, and B spectral reflectance to dominant wavelength based
  # on CIE chromaticity color space
  
  # see Wang et al 2015. MODIS-Based Radiometric Color Extraction and
  # Classification of Inland Water With the Forel-Ule
  # Scale: A Case Study of Lake Taihu
  
  require(colorscience)
  # chromaticity.diagram.color.fill()
  Xi <- 2.7689*R + 1.7517*G + 1.1302*B
  Yi <- 1.0000*R + 4.5907*G + 0.0601*B
  Zi <- 0.0565*G + 5.5943*B
  
  # calculate coordinates on chromaticity diagram
  x <-  Xi / (Xi + Yi +  Zi)
  y <-  Yi / (Xi + Yi +  Zi)
  z <-  Zi / (Xi + Yi +  Zi)
  
  # calculate hue angle
  alpha <- atan2( (x - 0.33), (y - 0.33)) * 180/pi
  
  # make look up table for hue angle to wavelength conversion
  cie <- cccie31 %>%
    mutate(a = atan2( (x - 0.33), (y - 0.33)) * 180/pi) %>%
    dplyr::filter(wlnm <= 700) %>%
    dplyr::filter(wlnm >=380)
  
  # find nearest dominant wavelength to hue angle
  wl <- cie[as.vector(sapply(alpha,function(x) which.min(abs(x - cie$a)))), 'wlnm']
  
  #out <- cbind(as.data.frame(alpha), as.data.frame(wl))
  
  return(wl)
}


## function for sampling either full stack or number if number is less than full stack.
sample_vals <- function (tbl, size, replace = FALSE, weight = NULL){
  ## assert_that(is.numeric(size), length(size) == 1, size >= 0)
  weight <- substitute(weight)
  index <- attr(tbl, "indices")
  sizes <- sapply(index, function(z) min(length(z), size)) # here's my contribution
  sampled <- lapply(1:length(index),
                    function(i) dplyr:::sample_group(index[[i]],  frac = FALSE, 
                                                     size = sizes[i],
                                                     replace = replace,
                                                     weight = weight))
  idx <- unlist(sampled) ## + 1
  grouped_df(tbl[idx, , drop = FALSE], vars = groups(tbl))
}

##############Bootstrapping functions

#Generate the bootsrap function
boot.med <- function(data, indices){
  dt<-data[indices,] %>% mutate(dummy = 'dummy')
  
  meds <- dt %>% group_by(dummy) %>% 
    summarise_at(vars(`1984`:`2018`), mean, na.rm = T) %>%
    select(-dummy)
  c(as.numeric(paste(meds[1,])))
}

## Create a follow up function to pull the mean and se (sd) of bootstrap iterations
boot.summary<- function(boots){
  summary <- tibble(mean = colMeans(boots$t),
                    se = apply(boots$t,2,sd),
                    bias = colMeans(boots$t) - boots$t0,
                    year = c(1984:2018))
  return(summary)
}


## Make look up table for AOD values
aodLUT <- function(path){
  df <- read.csv(path, stringsAsFactors = F) %>%
    mutate(COMID = as.character(COMID),
           year = year(date)) %>%
    filter(!is.na(blue)) %>%
    inner_join(lake.join %>% mutate(COMID = as.character(COMID)), by = c('COMID')) %>%
    distinct(system.index, date, .keep_all = T)
  
  lut <- df %>%
    select(COMID, date, lat = pour_lat, long = pour_long) %>%
    rowwise() %>%
    mutate(AOD = ifelse(ymd_hms(date) > ymd('2019/07/30'),
                        NA, pullMerra(date = date, lat = lat, long = long))) %>%
    ungroup()
  return(lut)
}
