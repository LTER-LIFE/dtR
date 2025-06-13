# ==============================================================================
# ==============================================================================
# Read data from KNMI, RWS, EMODnet bathymetry
# ==============================================================================
# ==============================================================================

## ==========================================
## weather (KNMI)
## ==========================================

read_KNMI <- function(file, dir="", attr=NULL, ...){
  KNMI <- NULL
  for (fi in file){
    fn <- paste(dir, fi, sep="/")
    fn <- gsub("//", replacement = "/", fn)
    TT <- readLines(fn)
    close(file(fn))

    skip <- grep("# STN", TT)-1                        # here the data start
    kk   <- read.csv(textConnection(TT), skip = skip[length(skip)], ...)
    
    CheckIt <- function(x){
      if (length(x) == 0) 
        x <- NA
      return(x)
    }
    
    WindSpeed      <- CheckIt(kk$FH/10      ) # m/s
    WindDirection  <- CheckIt(kk$DD         ) # degrees
    AirTemperature <- CheckIt(kk$T/10       ) # dgC
    Dewpoint       <- CheckIt(kk$TD/10      ) # dgC
    SolarRadiation <- CheckIt(kk$Q*1e4/3600 ) # J/m2/s (W/m2) was: J/cm2/hr
    
    if (length(SolarRadiation)==0) SolarRadiation <- NA
    Pressure       <- CheckIt(kk$P/10  )      # mbar (hPa) -> divide by 1000 for atm
    AirHumidity    <- CheckIt(kk$U/100  )      # -
    Cloud          <- CheckIt(kk$N     )
    SeaTemp        <- CheckIt(kk$TS/10 )
    
    Date <- as.POSIXct(strptime(paste(kk$YYYYMMDD, kk$HH, sep=" "),
                                format="%Y%m%d %H"))         # Date format
    
    WD             <- data.frame(station = kk[,1], 
                              datetime = Date, 
                              windspeed = WindSpeed, 
                              winddirection = WindDirection,
                              temperature = AirTemperature,
                              dewpoint = Dewpoint, radiation = SolarRadiation, 
                              pressure = Pressure, 
                              cloudcover=Cloud,
                              humidity = AirHumidity,
                              seatemperature = SeaTemp)
    
    KNMI <-  rbind(KNMI,   
                   WD)  #[1:(24*365),]
  }
  
  KNMI <- merge(KNMIstations[,c("station", "longitude", "latitude")],KNMI)
  vars <- data.frame(variable =c("windspeed", "winddirection",
                                 "temperature", "dewpoint", 
                                 "radiation", "pressure", 
                                 "cloudcover",
                                 "humidity", "seatemperature"),
                     description=c("Mean wind speed for the hour preceding the observation time stamp",
                                   "Mean wind direction for 10-minutes before time (360=N, 90=E, 0=calm, 990=variable)",
                                   "Temperature at 1.50 m at the time of observation",
                                   "Dew point temperature at 1.50 m at the time of observation",
                                   "Global radiation during the hourly division",
                                   "Air pressure reduced to mean sea level, at the time of observation",
                                   "Cloud cover (in octants), at the time of observation (9=sky invisible)",
                                   "Relative atmospheric humidity at 1.50 m at the time of observation",
                                   "Sea surface temperature at the time of observation"),
                     unit = c("m/s", "dg", "dgC", "dgC", "J/m2/s", "mbar (hPa)", "octants", "-", "dgC"),
                     original=c("FH", "DD", "T", "TD", "Q", "P", "N", "U", "TZ"))
  stats <- merge(KNMIstations, list(station = unique(KNMI$station)))
  ID    <- stats[,c("station", "longitude", "latitude")]
  names(ID)[1] <- "ID"
  
  attributes(KNMI)$variables  <- vars
  attributes(KNMI)$stations   <- stats
  attributes(KNMI)$ID         <- ID
  
  attributes(KNMI)$datasource <- "KNMI"
  attributes(KNMI)$file       <- file
  attributes(KNMI)$processing <-  paste("Created at", Sys.time())
  attributes(KNMI)$fun        <- "read_KNMI" 
  attributes(KNMI)$format     <- "wide" 
  if (length(attr))
    attributes(KNMI) <- c(attributes(KNMI), attr)
  
  class(KNMI) <- c("dtLife", class(KNMI))
  KNMI
}

## ==========================================
## water data (RWS)
## ==========================================

read_RWS <- function(file, dir="", attr=NULL, format="wide", ...){
  
  RWS <- NULL
  for (fi in file){
    fn <- paste(dir, fi, sep="/")
    fn <- gsub("//", replacement = "/", fn)
    RWS <- rbind(RWS, read.csv2(fn))  
  }
  RWS <- RWS[! is.na(RWS$X),]
  RWS <- RWS[! is.na(RWS$WAARNEMINGDATUM), ]
  select <- c("MEETPUNT_IDENTIFICATIE", "LOCATIE_CODE", 
              "WAARDEBEPALINGSMETHODE_OMSCHRIJVING",
              "WAARNEMINGDATUM", "WAARNEMINGTIJD..MET.CET.",  
              "PARAMETER_OMSCHRIJVING","GROOTHEID_OMSCHRIJVING", 
              "HOEDANIGHEID_OMSCHRIJVING",
              "EENHEID_CODE", "NUMERIEKEWAARDE", 
              "MEETAPPARAAT_OMSCHRIJVING", "X" ,      "Y" )

  cn <- c("stationName", "station", 
          "method",
          "date", "time",  "parameter", "description", 
          "property",  "unit", "value", 
          "sensor", "X" ,      "Y" )
  RWSdat <- RWS[, select]
  colnames(RWSdat) <- cn
  RWSdat <- RWSdat[RWSdat$value < 1e10, ]
  RWSdat$variable <- NA
  Vars <- rbind(
    c("Verzadigingsgraad",         "", "O2_percent", "Oxygen saturation percent"),   # NEEDS TO BE FIRST
    c("Biochemisch zuurstofverbruik met allylthioureum", "BZV5a", "BOD", "Biochemical oxygen demand"),
    c("chlorofyl-a",               "",      "Chl",     "Chlorophyll"),
    c("Chemisch zuurstofverbruik", "CZV",   "COD",     "Chemical oxygen demand"),
    c("Doorzicht",                 "",      "Secchi",  "Extinction depth"),
    c("Extinctie",                 "",      "Ext",     "Extinction- dimensionless"),
    c("fosfor totaal",             "Ptot",  "Ptot",    "Total phosphor"),
    c("koolstof",                  "Ctot",  "Ctot",    "Total carbon"),
    c("koolstof anorganisch",      "Cinorg", "PIC",    "Particulate inorganic carbon"),
    c("koolstof organisch",        "Corg",   "POC",    "Particulate organic carbon"),
    c("nitraat",                   "NO3",    "NO3",    "Nitrate"),  # PARAMETER_OMSCHRIJVING
    c("nitriet",                   "NO2",    "NO2",    "Nitrite"),
    c("ammonium",                  "NH4",    "NH4",    "Ammonium"),
    c("Onopgeloste stoffen",       "OS",     "SPM",    "Suspended particulate matter"),
    c("orthofosfaat",              "PO4",    "PO4",    "Phosphate"),
    c("siliciumdioxide",           "SiO2",   "SiO2",   "Silicium dioxide"),
    c("som nitraat en nitriet",    "sNO3NO2", "NOx",   "Nitrate+nitrite"),
    c("stikstof Kjeldahl",         "NKj",    "KjN",    "Kjehldahl nitrogen"),
    c("stikstof totaal",           "Ntot",   "Ntot",   "Total nitrogen"),
    c("Zuurgraad",                 "pH",     "pH",     "pH"),   # GROOTHEID_OMSCHRIJVING
    c("zuurstof",                  "",       "O2",     "Oxygen"),   
    c("sulfaat",                   "SO4",    "SO4",    "Sulphate"),
    c("waterstofcarbonaat",        "HCO3",   "HCO3",   "Bicarbonate"),
    c("Temperatuur",               "T",      "T",   "Temperature"),          
    c("Waterhoogte",               "WATHTE", "Height", "Waterheight"), 
    c("Saliniteit",                "SALNTT", "S",      "Salinity"),
    c("Waterhoogte berekend",      "WATHTBRKD", "Height", "Waterheight_estimated"))
  
  # Chloride, Chlorophyl, Doorzicht, Extinctiecoefficient
  # Opgelost organisch koolstof, Particulair organisch koolstof, Percentage zuurstof,
  # Saliniteit, Stikstof
  Vars <- as.data.frame(Vars)
  
    names(Vars) <- c("dutch_name", "dutch_code", "variable", "description")
    convfac=c(BOD=1,        COD=1,       
              Ptot=1e3/30.97376,    
              NO3=1e3/14.0067, NO2=1e3/14.0067, NH4=1e3/14.0067, 
              SPM=1, PO4=1e3/30.97376, SiO2=1e3/28.0855,       
              NOx=1e3/14.0067,  KjN=1e3/14.0067, Ntot=1e3/14.0067,
              pH=1,  SO4=1e3/32.065, 
              HCO3=1e3/(1+12.0107+3*15.9994),  Temp=1,      
              Height=0.01,     height_est =0.01, 
              O2=15.9994*2, 
              Secchi=0.1)  # was: dm
    newunits=c(BOD=NA, COD=NA, Ptot="mmol/m3", 
               NO3="mmol/m3", NO2="mmol/m3", NH4="mmol/m3", 
               SPM=NA, PO4="mmol/m3",       
               SiO2="mmol/m3", NOx="mmol/m3", KjN="mmol/m3", Ntot="mmol/m3",  
               pH="-",  
               SO4="mmol/m3", HCO3="mmol/m3",  # mg/L ???
               Temp=NA,      
               Height="m",     height_est ="m", O2="mmol/m3", Secchi="m")
    NN <- merge(data.frame(variable=names(convfac), convfac=convfac), 
                data.frame(variable=names(newunits), newunits=newunits), all=TRUE)
    
    VARS <- merge(Vars, NN, all=TRUE)
    for(i in 1:nrow(VARS)){
     
      ii <- which(RWSdat$description  == VARS$dutch_name[i])
      if (! length(ii))
       ii <- which(RWSdat$parameter == VARS$dutch_name[i])
      if (length(ii)){
        pars <- unique(RWSdat$description[ii])
        if (length(pars) > 1) 
          ii <- ii[which(RWSdat$description[ii] == "(massa)Concentratie")]
      }       
      if (length(ii)){
       RWSdat$description[ii] <- VARS$description[i]
       RWSdat$variable[ii]     <- VARS$variable[i]
       if (! is.na(VARS$newunits[i])){
        RWSdat$value[ii]    <- RWSdat$value[ii]*VARS$convfac[i]
        RWSdat$unit[ii]     <- VARS$newunits[i]
       }
      }
    }
    # manually clean the salinity, which also contains chlorinity!
    ii <- which(RWSdat$description == "Salinity" & RWSdat$value > 100)
    if (length(ii)) {
      RWSdat$description[ii] <- "Chlorinity"
      RWSdat$variable[ii] <- "Cl"
    }
    
  dd <- as.POSIXct(paste(RWSdat$date[1], RWSdat$time[1]), 
                   format="%d-%m-%Y %H:%M:%S")
  
  RWSdat$variable   [is.na(RWSdat$variable)]     <- "unknown"
  RWSdat$sensor     [is.na(RWSdat$sensor)]   <- "unknown"
  RWSdat$method     [is.na(RWSdat$method)]   <- "unknown"
  RWSdat$property   [is.na(RWSdat$property)] <- "unknown"
  RWSdat$unit       [is.na(RWSdat$unit)]     <- "unknown"
  RWSdat$description[is.na(RWSdat$description)] <- "unknown"
  if (!is.na(dd))
    RWSdat$datetime <- as.POSIXct(paste(RWSdat$date, RWSdat$time), 
                                  format="%d-%m-%Y %H:%M:%S")
  else 
    RWSdat$datetime <- as.POSIXct(paste(RWSdat$date, RWSdat$time), 
                                  format="%d/%m/%Y %H:%M:%S")
  Tfrom <- aggregate(as.Date(RWSdat$datetime), 
                  by=as.list(RWSdat[,c("variable","description", "property", "method", "unit", "sensor")]),
                  FUN=min, na.rm=TRUE)
  names(Tfrom)[ncol(Tfrom)] <- "from"
  Tto <- aggregate(as.Date(RWSdat$datetime), 
                  by=as.list(RWSdat[,c("variable","description", "property", "method", "unit", "sensor")]),
                  FUN=max)
  names(Tto)[ncol(Tto)] <- "to"
  variables <- merge(Tfrom, Tto)
  variables <- data.frame(variables, vnr=1:nrow(variables))
  
  wgs84 <- xy_to_wgs84(RWSdat$X, RWSdat$Y)
  RWSdat <- data.frame(wgs84, RWSdat)
  if (format == "wide")
    RWSdat <- data.frame(RWSdat, vnr=NA)  
  else
    RWSdat <- merge(RWSdat, variables)   ### THIS TAKES A VERY LONG TIME FOR BIG DATASETS
  
  stats <- unique(RWSdat[,c("stationName", "station", "X", "Y", "longitude", "latitude")])
  ID    <- stats[,c("station", "longitude", "latitude")]
  names(ID)[1] <- "ID"
  
  RWSdat <- RWSdat[,c("station", "longitude", "latitude", "datetime", "variable", "value", "unit", "vnr")]
  
  if (format == "wide"){
    RWSdat <- reshape(RWSdat[ , c(1:6)], 
                         direction = "wide", 
                         idvar     = c("station", "latitude", "longitude", "datetime"),   
                         timevar   = "variable")
    
    cn <- colnames(RWSdat)
    cn <- gsub("value.", "", cn)
    colnames(RWSdat) <- cn  
  } else format <- "long"
  
  attributes(RWSdat)$variables   <- variables
  attributes(RWSdat)$stations    <- stats
  attributes(RWSdat)$ID          <- ID
  attributes(RWSdat)$datasource  <- "RWS"
  attributes(RWSdat)$EPSG <- unique(na.omit(RWS$EPSG))
  attributes(RWSdat)$file <- file
  attributes(RWSdat)$processing <-  paste("Created at", Sys.time())
  attributes(RWSdat)$fun  <- "read_RWS" 
  attributes(RWSdat)$format <- format 
  if (length(attr))
    attributes(RWSdat) <- c(attributes(RWSdat), attr)
  class(RWSdat) <- c("dtLife", "data.frame")
  RWSdat
}

## ==========================================
## bathymetry (EMODnet)
## ==========================================

read_bathymetry <- function(file, dir="", attr=NULL, lonlim=NULL, latlim=NULL, 
                   levels=c( -5, 0, 5, 10, 20), by=1, ...){
  fn <- paste(dir, file, sep="/")
  fn <- gsub("//", replacement = "/", fn)
  
  morf <-  ncdf4::nc_open(fn)
  
  # read the elevation variable
  elev      <- ncdf4::ncvar_get(morf, "elevation")  
  
  # latitude and longitude
  latitude  <- ncdf4::ncvar_get(morf, "latitude")  
  longitude <- ncdf4::ncvar_get(morf, "longitude")  
  
  ncdf4::nc_close(morf)
  
  M.x <- 1: length(longitude)
  M.y <- 1: length(latitude)
  
  if (! is.null(lonlim))
    M.x <- which(longitude >= min(lonlim) & 
                 longitude <= max(lonlim))
  
  if (! is.null(latlim))
    M.y <- which(latitude >= min(latlim) & 
                 latitude <= max(latlim))
  
  if (by > 1) {
    M.x <- seq(from=M.x[1], to=M.x[length(M.x)], by=by)
    M.y <- seq(from=M.y[1], to=M.y[length(M.y)], by=by)
  }
  if (length(M.x) == 0 | length(M.y) == 0)
    stop ("cannot create bathymetry - no data selected")
  
  vars <- list(variable = c("longitude", "latitude", "depth"),
               unit     = c("dgE", "dgN", "m"),
               note = c("EPSG:4326", "original data: depth below the lowest astronomical tide, estimated as - elevation"))
    
  Bat <-  list(longitude=longitude[M.x], latitude=latitude[M.y], depth=-elev[M.x, M.y], 
               datasource="http://www.emodnet-bathymetry.eu/", file=file, 
               variables = vars, processing = paste("Created at", Sys.time()),
               EPSG=4326) 
  Bat$asp <- 1/cos((mean(Bat$latitude) * pi)/180)

  # convert NaN to NA
  Bat$depth[is.nan(Bat$depth)] <- NA
  
  # Create contourLines - remove too small polygons
  if (length(M.x) > 10 & length(M.y) > 10){
    cL <- with(Bat, contourLines(x=longitude, y=latitude, z=depth, 
                                    levels = levels))
    ll <- unlist(lapply(cL, FUN=function(x)length(x$x)))
    contours <- cL[ll>10]
  } else contours <- NULL 
  Bat$contours <- contours
  if (! is.null(attr)) Bat <- c(Bat, attr)
  attributes(Bat)$fun  <- "read_bathymetry"
  class(Bat) <- c("dtBathymetry", class(Bat))
  Bat
} 
