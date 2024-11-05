## ======================================================
##
## Working with FRRF and LABstaf data
##
## ======================================================

# ======================================================
# ======================================================
# Read FRRF data (act2, LabSTAF)
# ======================================================
# ======================================================

readFRRF <- function(file, dir = "",  # filename and directory
                     specs = NULL,    # specifications to be added to the data
                     attr = NULL,     # to be added to the attributes
                     txt = "csv",     # how to read the file
                     ...){          # other pars passed to the reading function
  if (txt == "csv") {
    Fread <- read.csv
    sep = ","
  }
  else if (txt == "csv2") {
    Fread <- read.csv2
    sep=";"
  }
  else if (txt %in% c("txt", "delim")) {
    Fread <- read.delim
    sep = "\t"
  }
  else if (txt == "delim2") {
    Fread <- read.delim2
    sep = "\t"
  }
  else if (txt == "table") {
    Fread <- read.table
    sep = ""
  }
  
  # filename
  if (dir == "") Fn <- file else Fn <- paste(dir, file, sep = "/")
  
  # read the file
  lines <- readLines(Fn, warn = FALSE)
  close(file(Fn))
  
  # determine if act2 or labstaf 
  act2     <- grep("Act2", lines) 
  LabSTAF  <- grep("LabSTAF", lines) 
  
  if (length(act2))
    ff <- parseFRRF.act2   (lines, file, specs, attr, Fread=Fread, 
                            sep = sep, ...)
  
  else if (length(LabSTAF))
    ff <- parseFRRF.LabSTAF(lines, file, specs, attr, Fread=Fread, 
                            sep = sep, ...)
  else if (length(grep("Acq", lines)))
    ff <- parseFRRF.acq(lines, file, specs, attr, Fread=Fread, 
                            sep = sep, ...)
  else
    stop ("Can only read act2 or LabSTAF files")
  
  ff
}

# ======================================

parseFRRF.acq <- function(lines, fn, specs = NULL, attr = NULL, 
                           Fread = read.csv,  sep = ",", 
                           ...){
  
  skip  <- grep("Acq", lines) - 1
  
  LL <- lines[1:(skip-1)]
  pars <- LL[-which(LL == "")]
  lines <- lines[-(1:skip)]
  
  Dat <- Fread(textConnection(lines), header = TRUE, ...)
  Dat[,1]     <- fn
  colnames(Dat)[1] <- "file"  
  Dat$file    <- fn
  Dat$Date <- paste(Dat$Date,Dat$Time)
  Dat$Time <- NULL
  cn <- colnames(Dat)
  oldcolnames <- c("Date", "X.Chl.", "Fo.or.F.", "Fm.or.Fm.", "Fv.Fm.or.Fq..Fm.", 
                   "JPSII.x.qJ",  "JPSII.x.qP" , "JPSII.x.qL", "dbar")
  newcolnames <- c("date", "Chl",  "Fo", "Fm", "Fv/Fm", "JPSII_qJ",  
                   "JPSII_qP" , "JPSII_qL", "E")
  
  ij <- match(oldcolnames, cn)
  
  cn[ij] <- newcolnames
  colnames(Dat) <- cn 
  
  ij <- which (colnames(Dat) %in% c("X.1", "X.2", "X.3"))
  if (length(ij)) 
    Dat <- Dat[, -ij]
  
  if (length(specs)){
    ELL <- as.data.frame(matrix(nrow = nrow(Dat), ncol = length(specs), 
                                data = specs, byrow = TRUE))
    names(ELL) <- names(specs)
    Dat <- cbind(Dat, ELL)
  }
  
  if (length(attr))
    attr(Dat, which = names(attr)) <- attr
  
  attributes(Dat)$system   <- "act2"
  attributes(Dat)$file     <- fn
  attributes(Dat)$date     <- date 
  attributes(Dat)$settings <- pars

  attributes(Dat)$processing <-  paste("Created at", Sys.time())
  attributes(Dat)$standardized <- FALSE
  
  Dat
}

# ======================================

parseFRRF.act2 <- function(lines, fn, specs = NULL, attr = NULL, 
                           Fread = read.csv,  sep = ",", 
                           ...){

  getFRRFdate <- function(lines){
    idate <- grep("File date", lines)
    itime <- grep("File time", lines)
    date  <-  paste(strsplit(lines[idate], sep)[[1]][2], 
                    strsplit(lines[itime], sep)[[1]][2])
    lines <<- lines[-(1:(itime+1))]
    return(date)
  } 
  
  getFRRFfit <- function(lines){
    ip <- grep("Alpha:", lines)-1
    il <- grep("SErP:", lines)
    fits  <-  Fread(textConnection(lines[ip:il]), header = TRUE, ...)
    lines <<- lines[-((ip-1):(il+1))]
    return(fits)
  }    
  
  date  <- getFRRFdate(lines)
  fit   <- getFRRFfit(lines)
  skip  <- grep("Saq", lines) - 1
  
  LL <- lines[1:(skip-1)]
  LL <- LL[-which(LL == "")]
  pars  <- Fread(textConnection(LL), header = FALSE, ...)
  lines <- lines[-(1:skip)]
  
  Dat <- Fread(textConnection(lines), header = TRUE, ...)
  Dat[,1]     <- fn
  colnames(Dat)[1] <- "file"  
  Dat$date    <- date
  Dat$file    <- fn
  cn <- colnames(Dat)
  oldcolnames <- c("X.Chl.", "rP", "rP.1", "F.", "Fm.", "Fq..Fm.")
  newcolnames <- c("Chl", "rP_measured", "rP_fitted", "Fo", "Fm", "Fv/Fm")

  ij <- match(oldcolnames, cn)
  
  cn[ij] <- newcolnames
  colnames(Dat) <- cn 
  
  ij <- which (colnames(Dat) %in% c("X.1", "X.2", "X.3"))
  if (length(ij)) 
    Dat <- Dat[, -ij]
  
  if (length(specs)){
    ELL <- as.data.frame(matrix(nrow = nrow(Dat), ncol = length(specs), 
                                data = specs, byrow = TRUE))
    names(ELL) <- names(specs)
    Dat <- cbind(Dat, ELL)
  }
  
  if (length(attr))
    attr(Dat, which = names(attr)) <- attr
  
  attributes(Dat)$system <- "act2"
  attributes(Dat)$file   <- fn
  attributes(Dat)$date   <- date 
  attributes(Dat)$settings <- pars
  attributes(Dat)$fitted.pars <- fit 
  attributes(Dat)$processing <-  paste("Created at", Sys.time())
  attributes(Dat)$standardized <- FALSE
  
  Dat
}

# ======================================

parseFRRF.LabSTAF <- function(lines, fn, specs = NULL, attr = NULL, 
                              Fread = read.csv, sep = ",", ...){
  
  getaLHII <- function(lines){
    i1 <- grep("aLHII", lines)[1]
    strsplit(lines[i1], ":")[2]
  }
  getKa <- function(lines){
    i1 <- grep("Ka", lines)[1]
    strsplit(lines[i1], ":")[2]
  }  
  getFRRFsetup <- function(lines){
    i1 <- grep("STAF", lines)[1]
    i2 <- grep("Fo", lines)[1] - 1
    setup  <-  gsub(sep, " ", lines[i1:i2])
    setup <- gsub("<b5>", "micro ", setup)
    setup <- gsub("<b0>", "dg ", setup)
    setup <- gsub("<b2>", "^2 ", setup)
    lines <<- lines[-(1:i2)]
    return(setup)
  }    
  getFRRFdate <- function(lines){
    idate <- grep("Date", lines)
    itime <- grep("Time", lines)
    date  <-  paste(strsplit(lines[idate], sep)[[1]][2], 
                    strsplit(lines[itime], sep)[[1]][2])
    lines <<- lines[-c(idate,itime)]
    return(date)
  }    
  getFRRFfit <- function(lines){
    ip <- grep("Alpha:", lines)-1
    il <- grep("GOPIIm:", lines)
    if (! length(il)) il <- grep("rPm:", lines)
    fits  <-  gsub(sep, " ", lines[ip:il])
    fits <- gsub("<b5>", "micro ", fits)
    fits <- gsub("<b0>", "dg ", fits)
    fits <- gsub("<b2>", "^2 ", fits)
    fits <- gsub("<b3>", "^3 ", fits)
    lines <<- lines[-((ip-1):(il+1))]
    return(fits)
  }    
  aLHII_0 <- getaLHII(lines) 
  Ka     <- getKa(lines) 
  date   <- getFRRFdate(lines)
  fit    <- getFRRFfit(lines)
  pars   <- getFRRFsetup(lines)

  nend <- grep("Gaq", lines)[1]-1
  Dat  <- Fread(textConnection(lines[2:nend]), header = TRUE, ...)
  Dat[,1] <- fn
  colnames(Dat)[1] <- "file"  
  Dat$date    <- date
  Dat$file    <- fn
  cn <- colnames(Dat)
  oldcolnames <- c("F.", "Fm.")
  newcolnames <- c("Fo", "Fm")
  
  ij <- match(oldcolnames, cn)
  
  cn[ij] <- newcolnames
  colnames(Dat) <- cn 
  
  
  if (length(specs)){
    ELL <- as.data.frame(matrix(nrow = nrow(Dat), ncol = length(specs), 
                                data = specs, byrow = TRUE))
    names(ELL) <- names(specs)
    Dat <- cbind(Dat, ELL)
  }
  if (length(attr))
    attr(Dat, which=names(attr)) <- attr
  
  attributes(Dat)$system <- "LabSTAF"
  attributes(Dat)$file <- fn
  attributes(Dat)$date <- date 
  attributes(Dat)$aLHII_0 <- aLHII_0 
  attributes(Dat)$Ka      <- Ka
  attributes(Dat)$settings <- pars
  attributes(Dat)$fitted.pars <- fit 
  attributes(Dat)$processing <- paste("Created at", Sys.time())
  attributes(Dat)$standardized <- TRUE
  
  Dat
}

# ======================================================
# ======================================================
# Standardize FRRF data - using the absorption method
# ======================================================
# ======================================================

standardizeFRRF <- function(frrf, 
                            Fblanc = 0,   # background fluorescence in water 
                            Ka = 11800,   # instrument-specific constant, units of /m
                            convJVPII = 3600/1000,# from micromol/s to mmol/hour
                            aLHII_0 = NA, # light harvesting 
                            na.rm = TRUE, # to remove the failed data points or not
                            verbose = FALSE 
                            ){  
  
  
  # save the attributes
  ATT <- attributes(frrf)
  ATT <- ATT[!names(ATT) %in% c("names", "row.names", "class")]
  
  if (is.null(ATT$system))
    ATT$system <- "unknown"
  if (is.null(ATT$standardized))
    ATT$standardized <- FALSE
  
  if (ATT$system == "LabSTAF") {
    if (verbose)
       warning ("*RE* standardizing frrf of type LabSTAF !!!  this is already done") 
    if (is.null(frrf$Fo_uc)){
      frrf$Fo_uc    <- frrf$Fo
      frrf$Fm_uc    <- frrf$Fm
      frrf$JVPII_uc <- frrf$JVPII
    }
  }
  
  if (verbose & ATT$standardized) 
    warning ("frrf file is already standardized-will standardize again") 
  
  # add blanc
  ff         <- frrf
  ff$Fblanc  <- Fblanc
  
  # keep uncorrected values (uc)
  if (! ATT$standardized){
    ff$Fo_uc    <- ff$Fo
    ff$Fm_uc    <- ff$Fm
    ff$JVPII_uc <- ff$JVPII
  }
  
  # remove failed data points
  if (na.rm){
    ii <- which (is.na(ff$Fo))
    if (length(ii)) 
      ff <- ff[-ii, ]
  }
  
  if (nrow(ff) == 0) return(NULL)
  # correct for the blanc (fluorescence in the absence of Chl)
  ff$Fo      <- with(ff, Fo_uc - Fblanc) 
  ff$Fm      <- with(ff, Fm_uc - Fblanc)
  
  # calculate values using the absorption method
  ff$Fq      <- with(ff,  Fm - Fo )
  
  # cross-sectional surface of PSII system
  ff$a_LHII  <- with(ff,  Fo * Fm / Fq * (Ka / 10^6))     
  ff$FqFm    <- with(ff,  Fq/ Fm  )
  
  # volumetric e-flux, [convJPII * micromol m-3 s-1]
  # note: we need a_LHII at E=0.
  if (is.na(aLHII_0))
    aLHII_0 <- mean(ff$a_LHII[ff$E <= 0])
  
  Noa <- (length(aLHII_0) == 0)
  if (! Noa) Noa <- (is.na(Noa) | is.nan(aLHII_0) | is.infinite(aLHII_0))
    
  if (Noa) {# If no data at E=0; use linear regression
    ffS <- ff [ff$E < 100,]
    ffS <- ffS[! is.infinite(ffS$a_LHII), ]
    if (nrow(ffS) > 1)
      aLHII_0 <- coef(lm(ffS$a_LHII ~ ffS$E))[1]
    else if (nrow(ffS) == 1)
      aLHII_0 <- ffS$a_LHII
    else 
      aLHII_0 <- NA
   }   
  
  ff$JVPII   <- ff$Fq/ ff$Fm* aLHII_0 * ff$E * convJVPII 
  
  # save attributes
  attributes(ff) <- c(attributes(ff), ATT)
  
  attributes(ff)$processing <- c(attributes(ff)$processing, 
                paste("Standardized with Fblanc = ", unique(Fblanc), "at", Sys.time()),
                paste("JVPII calculated (absorption method), with conversion factor = ", 
                      convJVPII, "at", Sys.time()))
  attributes(ff)$standardized <- TRUE
  if (ATT$system == "FRRF") 
    attributes(ff)$check <- checkFRRF(ff)  
  if (convJVPII == 1) 
    U_JVPII <- "umol photons / m3/s"
  else if (convJVPII == 3.6) 
    U_JVPII <- "mmol photons/m3/hour"
  else if (convJVPII == 86.4) 
    U_JVPII <- "mmol photons/m3/day"
  else 
    U_JVPII <- paste("umol photons/m3/", convJVPII, "s", sep="")
  attributes(ff)$unit_JVPII <- U_JVPII
  
  attributes(ff)$aLHII_0    <- aLHII_0
  attributes(ff)$Ka         <- Ka
  return(ff)
}  

# ======================================================
# Check the data - flag suspicious values
# ======================================================

flagFRRF <- function(frrf, 
                     QA_ADC    = c(30, 80), 
                     QA_RSigma = c(0.035, 0.07),
                     QA_QR     = c(6, NA))
{
  frrf$flag <- ""
  iADC <- which(frrf$ADC < QA_ADC[1] |  
                frrf$ADC > QA_ADC[2])
  frrf$flag[iADC] <- paste(frrf$flag[iADC],"A",sep="")


  iR <- which(frrf$RSigma < QA_RSigma[1] | 
              frrf$RSigma > QA_RSigma[2])
  frrf$flag[iR] <- paste(frrf$flag[iR],"R",sep="")
  
  iQ <- which(frrf$QR < QA_QR[1])
  frrf$flag[iQ] <- paste(frrf$flag[iQ],"Q",sep="")

  attributes(frrf)$processing <- c(attributes(frrf)$processing, 
                                   paste("flagged at", Sys.time()))
  
  frrf
}

# ======================================================

checkFRRF <- function(frrf, 
                      QA_ADC    = c(30, 80), 
                      QA_RSigma = c(0.035, 0.07),
                      QA_QR     = c(6, NA),
                      plotit    = FALSE)
{
  check <- data.frame()
  frrf$flag <- ""
  iADC <- which(frrf$ADC < QA_ADC[1] |  
                frrf$ADC > QA_ADC[2])
  out <- length(iADC)
  
  check <- data.frame(name = "ADC", OK = out==0,
                      minValue = min(frrf$ADC), 
                      maxValue = max(frrf$ADC),
                      minQA= QA_ADC[1], 
                      maxQA=QA_ADC[2], 
                      exceed = out)
  
  
  iR <- which(frrf$RSigma < QA_RSigma[1] | 
              frrf$RSigma > QA_RSigma[2])
  frrf$flag[iR] <- paste(frrf$flag[iR],"R",sep="")
  out <- length(iR)
  check <- rbind(check, data.frame(name = "Rsigma", 
                                   OK = out==0,
                                   minValue = min(frrf$RSigma), 
                                   maxValue = max(frrf$RSigma),
                                   minQA= QA_RSigma[1], 
                                   maxQA=QA_RSigma[2], exceed = out))
  
  iQ <- which(frrf$QR < QA_QR[1])
  
  frrf$flag[iQ] <- paste(frrf$flag[iQ],"Q",sep="")
  out <- length(iQ)
  check <- rbind(check, data.frame(name = "QR", 
                                   OK = out==0,
                                   minValue = min(frrf$QR), 
                                   maxValue = max(frrf$QR),
                                   minQA= QA_QR[1], 
                                   maxQA=NA, exceed = out))
  
  check
}

