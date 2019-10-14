## Function for dividing data up 60/40 into train/test data

traintest <- function(x, ungroup = T){
  
  # thresh <- ifelse(grepl('doc', x$siteParam[1]) == T, 20, 
  #                  ifelse(grepl('chl', x$siteParam[1]) == T, 30,
  #                         ifelse(grepl('secchi', x$siteParam[1]) == T, 5, 60)))
  # x <- x %>%
  #   mutate(event = ifelse(value >= thresh, 1, 0))
  
  ##  We also find a mean value for each lake and bin into quantiles and make sure each
  ## quantile is represented in each fold to make them representative. Basically takes a 
  ## more stratified sample
  if(ungroup == F){
    quantiles <- x %>%
      group_by(COMID) %>%
      summarise(mean.v = mean(value)) %>%
      mutate(mag = cut(mean.v, quantile(x = mean.v,
                                        seq(0,1,0.1)), include.lowest = T),
             mag = factor(mag, labels = seq(.1,1,.1)))
    
    train <- x %>%
      left_join(quantiles, by = 'COMID')
    
    set.seed(22)
    out <- partition(train, p = .6, id_col = 'COMID', cat_col = 'mag', list_out = F)
  }
  ## To just stratify by observed value.  Smaller increments at higher quantiles ensure extreme events appear
  ## in both train and test data
  if(ungroup == T){
    x <- x %>%
      mutate(mag = cut(value, quantile(x = value,
                                       c(0, 0.2, 0.4, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.93, 0.96,1),
                                       include.lowest = T)),
             mag = factor(mag, 
                          labels = c(0.2, 0.4, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.93, 0.96, 1)))

    set.seed(22)
    train <- x %>%
      group_by(mag, sat) %>%
      sample_frac(.6) %>%
      ungroup() %>%
      mutate(.partitions = 1)
    
    test <- x %>%
      filter(!(uniqueID %in% train$uniqueID)) %>%
      mutate(.partitions = 2)
    
    out <- train %>%
      bind_rows(test)
  }
  return(out)
}

#Function for pulling out different groups of input features
featSelect <- function(numAdd, param, vars){
  if(numAdd == 0){
    features <- c('blue','green','red', 'nir','swir1', 'swir2', "NR", "BR", "GR",
                  "SR", "BG", "BN", "BS", 'GS', 'GN', 'fai', 'ndvi', 'ndwi', 
                  'sat', 'hillshade', 'zenith')
  }else if(numAdd > 0){
    features = vars %>%
      filter(parameter == param) %>%
      select(var) %>%
      str_split(pattern = ', ') %>%
      first()
    
    features = c(features[!features %in% c('blue','green','red', 'nir','swir1',
                                           'swir2', "NR", "BR", "GR", "SR", "BG",
                                           "BN", "BS", 'GS', 'GN', 'fai', 'ndvi',
                                           'ndwi', 'sat', 'hillshade', 'zenith')][1:paste0(numAdd)],
                 'blue','green','red', 'nir','swir1', 'swir2', "NR", "BR", "GR",
                 "SR", "BG", "BN", "BS", 'GS', 'GN', 'fai', 'ndvi', 'ndwi',
                 'sat', 'hillshade', 'zenith')
  }else if(numAdd == -99){
    features = vars %>%
      filter(parameter == param) %>%
      select(withinSD) %>%
      str_split(pattern = ', ') %>%
      first()
    #New method taking just the first 3 non-optical features from the RF RFE.
    features = c(features,'blue','green','red', 'nir','swir1', 'swir2', "NR", "BR", "GR",
                 "SR", "BG", "BN", "BS", 'GS', 'GN', 'fai', 'ndvi', 'ndwi',
                 'sat', 'hillshade', 'zenith')
  }else if(numAdd == -9999){
    features = vars %>%
      filter(parameter == param) %>%
      select(withinSD2) %>%
      str_split(pattern = ', ') %>%
      first()
    #New method taking just the first 3 non-optical features from the RF RFE.
    features = c(features[features!= 'PctWdWet'],'blue','green','red', 'nir','swir1', 'swir2', "NR", "BR", "GR",
                 "SR", "BG", "BN", "BS", 'GS', 'GN', 'fai', 'ndvi', 'ndwi',
                 'sat', 'hillshade', 'zenith')
  }
  return(features)
}



## Function for building XgBoost Models
boost <- function(param, data, vars, numAdd,  filter = T, weight = F, log = 'none', save = T){
  
  data.flat <- data %>%
    filter(parameter == param) %>%
    select(tt) %>%
    flatten_df()

  if(filter == T){
    data.flat <- data.flat %>%
      filter(dswe_sd == 0)
  }
  
  if(log == 'ln'){
    data.flat <- data.flat %>%
      mutate(value = log(value)) %>%
      filter(is.finite(value))
  }
  
  if(weight == T){
    data.flat <- data.flat %>%
      mutate(weights = (value - min(value))/(max(value)-min(value)),
             weights = ifelse(value < quantile(value, .6), 1, 2.5))
  }
  
  train <- data.flat %>%
    filter(.partitions == 1) %>%
    select(-.partitions)
  
  test <- data.flat %>%
    filter(.partitions == 2) %>%
    select(-.partitions)
  
  features <- featSelect(numAdd, param, vars)

  # One hot encoding for categorical features
  encoder <- train %>% select(features) %>% 
    onehot(.,max_levels = 12) 
  
  train.dat <- predict(encoder, train %>%
                         select(features))
  
  train.label <- train$value
  test.label <- test$value

  train_Data <- xgb.DMatrix(data = train.dat, label = train.label)
  
  if(weight == T){
    train_Data <- xgb.DMatrix(data = train.dat, label = train.label, weight = train$weights)
  }

  test.dat <- predict(encoder, test %>% select(features))
  
  test_Data <- xgb.DMatrix(data = test.dat, label = test.label)
  
  #Set Conservative hyperparametrs
  xgb.parameters<-list(
    objective =  "reg:linear",
    eval_metric = 'rmse',
    booster = "gbtree",
    max_depth = 3,
    eta = .2,
    gamma = .1, 
    subsample = .4,
    colsample_bytree = .4,
    min_child_weight = 1,
    lambda = 3,
    max_delta_step = 0 #This is supposed to help with eneven data distribution (e.g. very few high obs)
  )
  
  Training <- xgb.train(params = xgb.parameters,
                        data = train_Data,
                        nrounds = 300,
                        watchlist = list(train = train_Data, test = test_Data),
                        verbose = TRUE,
                        print_every_n = 20,
                        early_stopping_rounds = 20,
                        maximize = F,
                        missing = -9999)
  
  #This is just to look at the hueristics of training, come back to it later.
  eval.log <- Training$evaluation_log %>%
    mutate(Param = param, Method = 'xgboost', log = log)
  
  if(log == 'ln'){
    pred <- exp(predict(Training, test_Data))
    actual <- exp(test.label)
  }else{
    pred <- predict(Training, test_Data)
    actual <- test.label
  }

  #err <- as.numeric(sum(as.integer(pred > 0.5) != test.label))/length(test.label)
  #print(paste("test-error=", err))
  
  uniqueID <- test$uniqueID
  SiteID <- test$SiteID
  
  feat.names <- dimnames(train.dat)[[2]]
  feat.imp <- xgb.importance(feat.names, model = Training)
  
  if(save == T){
  saveKey <- paste0(iteration,'_',log,'_', param)
  xgb.save(Training, paste0('models/',saveKey,'.model'))
  feat.imp <- feat.imp %>%
    mutate(Parameter = param,
           log = log)
  write_feather(feat.imp, paste0('out/featureImp/',saveKey,'.feather'))
  }
  
  results <- tibble(Param = param, Method = 'xgboost', 
                    Predicted = pred, Actual = actual, log = log, 
                    uniqueID = uniqueID, SiteID = SiteID)
  
  ###
  #confusionMatrix(factor(results$Predicted), factor(results$Actual))
  
  return(list(results = results, eval.log = eval.log))
}


## Function for mapping over evaluation lakes data and generating predictions
EvalPreds <- function(id, paths, sites.eval, param, log, eval = T, version = 'new', model = NULL, features = NULL){
  
  if(version == 'old'){
    model = xgb.load(paste0('../aquaModel/models/reg_StrictSecchiSD_UnWeightedln_none_',param,'.model'))
  }else if(is.null(model)){
    saveKey <- paste0(iteration,'_',log,'_', param)
    model <- xgb.load(paste0('models/',saveKey,'.model'))}
  
  if(is.null(features)){
  inputs <- featSelect(numAdd, param, inputs)}

  
  # Really buggy but no serious ones, set it up with tryCatch so it'll return an empty data frame if it fails
  predicted <- tibble(system.index = NA, COMID = NA, date = NA, year = NA, month = NA, value = NA, parameter = NA)
  try({
    df <- read.csv(grep(paste0('/',id,'.csv'),paths, value = T), stringsAsFactors = F) %>%
      select_at(vars(-contains('hylak_id'))) %>%
      mutate(COMID = as.character(COMID),
             year = year(date),
             period = ifelse(year < 2009, 1, 2)) %>%
      filter(!is.na(blue)) %>% #,
             #hillshadow == 1) %>%
      inner_join(sites.eval %>% mutate(COMID = as.character(COMID)), by = c('COMID', 'period')) %>%
      distinct(.keep_all = T)
    
    pix.cut <- max(df$pixelCount)/2
    
    df <- df %>%
      mutate(nir = ifelse(sat == '8', 1.01621*nir + 99.28166, nir),
             red = ifelse(sat == '8', 0.9346232*red + 83.54783, red),
             green = ifelse(sat == '8', 0.8975912*green + 101.77444, green),
             NR = nir/red,
             clouds = Cloud_Cover,
             BR = blue/red,
             GR = green/red,
             NR = nir/red,
             SR = swir1/red,
             BG = blue/green,
             BN = blue/nir,
             BS = blue/swir1,
             GS = green/swir1,
             GN = green/nir,
             fai = nir - (red + (swir1-red)*((830-660)/(1650-660))),
             ndvi = ((nir-red)/(nir+red)),
             ndwi = ((green- swir1)/(green + swir1)),
             sat = factor(sat, levels = c('5','7','8')),
             #date = ymd_hms(date),
             month = as.numeric(month(date)),
             fui.hue = fui.hue(red, blue, green),
             #brightness = rgb2hsv(r=red,g= green, b= blue, maxColorValue = 10000)[3],
             nir = ifelse(sat == '8', nir + 100, nir)) %>%
      filter(pixelCount > pix.cut,
             pixelCount > 8)


    # One hot encoding for categorical features
    if(version == 'old'){
      df <- df %>% mutate(sat = factor(sat, levels = c('5','7','8'), labels = c('LT05', 'LE07', 'LC08')))
    }
    
    if('AOD' %in% features){
      lut <- read_feather('out/aodLUT_NLA2012.feather')
      
      df <- df %>%
        left_join(lut)
    }
    
    if(is.null(features)){
    encoder <- df %>% select(inputs) %>% 
      onehot()
    lake.input <- predict(encoder, df %>% select(inputs))
    }else if(eval == T){encoder <- df %>% select(-c(COMID, system.index, date, siteParam, parameter, latR, longR)) %>% onehot()
    lake.input <- predict(encoder, df) %>% as.data.frame() %>% select(features) %>% as.matrix()
    }else if(eval == F & "OBJECTID" %in% names(df)){
      df <- df %>% rename_at(vars(OBJECTID:ONOFFNET), tolower)
      facts <- c(names(df %>% select_if(is.factor)), '.geo', 'fdate', 'date', 'system.index', 'COMID')
      facts <- facts[!facts %in% c('sat', 'inStreamCat')]
      encoder <- df %>% select(-facts)
      encoder <- onehot(encoder)
      lake.input <- predict(encoder, df) %>% as.data.frame() %>% select(features) %>% as.matrix()
    }else{encoder <- df %>% 
      select(-c(COMID, system.index, date, grand_id, .geo, inputs:gnis_name, reachcode:poly_src)) %>%
      mutate(lake_type = factor(lake_type)) %>% onehot()
    lake.input <- predict(encoder, df) %>% as.data.frame() %>% select(features) %>% as.matrix()}
    
    
    if(log == 'ln'){
      value = tibble(value = exp(predict(model, lake.input)))
    }else{
      value = tibble(value = predict(model, lake.input))
    }
    
    if(eval == T){
      predicted <- df %>%
        select(system.index, COMID, date, year, month, paramField = parameter) %>%
        bind_cols(value) %>%
        mutate(parameter = param) %>%
        filter(parameter == paramField) %>%
        select(-paramField)}
    if(eval == F){
      predicted <- df %>%
        select(system.index, COMID, date, year, month, region) %>%
        bind_cols(value) %>%
        mutate(parameter = param)
    }
  })
  return(predicted)
}

#Function for making figures over space and time
spaceTimeFigs <- function(con, abs = T){
  if(abs == T){
    ggplot(errorSum %>% filter(Metric == 'mae'|Metric == 'bias', Param == con), aes(x = fct_reorder(quantLabs, order), y = Error, color = Metric, group = interaction(Metric, Variable, log), linetype = log )) + 
      geom_point() + 
      geom_line() + 
      facet_wrap(~Variable, scales = 'free')  +  
      theme_bw() +
      #scale_y_continuous(labels = scales::percent) +
      theme(axis.text.x = element_text(angle = 90, vjust = .5),
            legend.position = 'bottom') +
      labs(x = 'Quantile', y = 'Error (m)', title = paste0(con,' Absolute Errors'))
  }else{
    ggplot(errorSum %>% filter(Metric == 'smape'|Metric == 'p.bias', Param == con), aes(x = fct_reorder(quantLabs, order), y = Error, color = Metric, group = interaction(Metric, Variable, log), linetype = log)) + 
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


## Write a function to sample lakes and run the stl decomp.
bootstrap.stl <- function(iteration, full.preds, ice.months, sample.size = 100){
  
  ## Create dummy variable
  dummy <- expand.grid(year = seq(1985,2018), 
                       month = c(1:12), 
                       region = unique(region$region),
                       parameter = c('tss', 'chl_a', 'secchi')) %>% 
    mutate(iceID = paste0(month,region)) %>%
    right_join(ice.months %>% 
                 select(month,region, num.vis) %>% 
                 mutate(iceID = paste0(month,region)))
  
  ## Sample lakes in each HUC
  if(is.null(sample.size)){
    samp <- full.preds %>%
      distinct(COMID, .keep_all = T) %>%
      group_by(region) %>%
      nest() %>%
      mutate(n = purrr::map(data, nrow),
             samp = purrr::map2(data, n, sample_n, replace = T)) %>%
      unnest(samp)
  }else{samp <- full.preds %>%
    distinct(COMID, .keep_all = T) %>%
    group_by(region) %>%
    sample_n(sample.size, replace = T) %>%
    select(COMID, region)}
  
  ## Calculate median based on those lakes  
  yearly.median <- full.preds %>%
    na.omit() %>%
    filter(COMID %in% samp$COMID,
           year < 2019) %>%
    group_by(COMID, year, month, region, parameter) %>%
    dplyr::summarize(value = median(value, na.rm = T)) %>%
    group_by(year, month, region, parameter) %>%
    dplyr::summarize(median = median(value, na.rm = T),
                     sd = sd(value, na.rm =T)) %>%
    right_join(dummy)%>%
    arrange(parameter, region, year, month)%>%
    filter(parameter == 'secchi') %>% # Remove if adding more params than secchi
    mutate(date = as.POSIXct(paste0(year,'/',month,'/',01)))
  
  ## Run stl and save results in a data frame
  
  stl <- yearly.median %>%
    group_by(region, parameter) %>%
    nest() %>%
    mutate(stl = future_map(data, ~stlplus(x = .$median, t = .$date, 
                                           n.p = .$num.vis[1], s.window = 25,#OG 25
                                           sub.labels = c(1:.$num.vis[1]))),
           date = purrr::map(stl, 'time'),
           stl = purrr::map(stl, 'data')) %>%
    dplyr::select(-data) %>%
    unnest() %>%
    mutate(year = year(date),
           month = month(date)) %>%
    left_join(region %>% st_set_geometry(NULL)) %>%
    mutate(iteration = iteration)
  
  return(stl)
} 

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




# We want ~ two hundred lakes per huc, but a couple have fewer than 200.  This function take
# either all the lakes if less than sample size or sample size.
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
