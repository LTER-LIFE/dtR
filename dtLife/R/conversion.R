
# -------------------------------------------------------------------------
# spatial conversion and manipulation functions
# -------------------------------------------------------------------------

xy_to_wgs84 <- function(X, Y, EPSG=25831){
  
  # WGS84=EPSG4326
  SS <- data.frame(X=X, Y=Y)
  str <- paste("EPSG:", EPSG, sep="") 
  try(
    DD <- st_as_sf(SS, coords=c("X","Y"), crs=str), 
    silent=TRUE)
  dt <- st_coordinates(st_transform(DD, "EPSG:4326"))
  
  colnames(dt) <- c("longitude", "latitude")
  dt
}

# -------------------------------------------------------------------------

wgs84_to_xy <- function(longitude, latitude, EPSG=25831){
  
  # WGS84=EPSG4326
  SS <- data.frame(longitude=longitude, latitude=latitude)
  try(
    DD <- st_as_sf(SS, coords=c("longitude","latitude"), crs="EPSG:4326"), 
    silent=TRUE)
  str <- paste("EPSG:", EPSG, sep="") 
  dt  <- st_transform(DD, str)
  st_coordinates(dt)
}

# -------------------------------------------------------------------------
#FROM: F.H. Schreutelkamp, & ir. G.L. Strang van Hees,
# Benaderingsformules voor de transformatie
#tussen RD- en WGS84-kaartcoördinaten

#' # Martinitoren Groningen:
#' rd_to_wgs84(233883.131, 582065.167)

rd_to_wgs84 <- function(X, Y) {
  X0       <- 155000  # base coordinates of RD
  Y0       <- 463000
  
  Coeff <- list(A0=663304.11, B0=5780984.54, 
                A1=99947.539, B1=3290.106, 
                A2=20.008,    B2=1.310, 
                A3=2.041,     B3=0.203, 
                A4=0.001,     B4=0.000)
  # for zone 32
  Coeff2 <- list(A0=252878.65, B0=5784453.44, 
                A1=99919.783,  B1=-4982.166, 
                A2=-30.208,    B2=3.016, 
                A3=2.035,      B3=-0.309, 
                A4=-0.002,     B4=0.001)
  
  dX <- (X - X0) * 10^-5
  dY <- (Y - Y0) * 10^-5
  dX2 <- dX*dX
  dX3 <- dX2*dX
  dX4 <- dX2*dX2
  
  dY2 <- dY*dY
  dY3 <- dY2*dY
  dY4 <- dY2*dY2
  
  longitude <- with (Coeff,
             A0 + A1*dX - B1*dY + A2*(dX2 - dY2) - B2*(2*dX*dY) +
               A3*(dX3 - 3*dX*dY2) - B3*(3*dX2* dY - dY3) +
               A4*(dX4 - 6*dX2*dY2 + dY4) - B4*(4*dX3*dY - 4*dY3*dX))
  
  latitude <- with (Coeff, 
             B0 + B1*dX + A1*dY + B2*(dX2 - dY2) + A2*(2*dX*dY) +
               B3*(dX3 - 3*dX*dY2) + A3*(3*dX2*dY - dY3) +
               B4*(dX4 - 6*dX2*dY2 + dY4) + A4*(4*dX3*dY - 4*dY3*dX))
  data.frame(longitude=longitude*1e-5, latitude=latitude*1e-5)
}

# -------------------------------------------------------------------------

wgs84_to_rd <- function(longitude, latitude) {
  E0   <- 663304.11e-5
  N0   <- 5780984.54e-5

  Coeff <- list(C0=155000,    D0=463000,
                C1=99944.187, D1=-3289.996, 
                C2=-20.039,   D2=0.668,
                C3=-2.042,    D3=0.066, 
                C4=0.001,     D4=0.000)
  Coeff2 <- list(C0=155000,     D0=463000,
                 C1= 99832.079, D1=4977.793, 
                 C2=30.280,     D2=1.514,
                 C3=-2.034,     D3=-0.099, 
                 C4=-0.001,     D4=0.000)
 
  dE <- (longitude - E0) 
  dN <- (latitude - N0) 
  dE2 <- dE*dE
  dE3 <- dE2*dE
  dE4 <- dE2*dE2
  
  dN2 <- dN*dN
  dN3 <- dN2*dN
  dN4 <- dN2*dN2
  
  X <- with (Coeff,
             C0 + C1*dE - D1*dN + C2*(dE2 - dN2) - D2*(2*dE*dN) +
               C3*(dE3 - 3*dE*dN2) - D3*(3*dE2*dN - dN3) +
               C4*(dE4 - 6*dE2*dN2 + dN4) - D4*(4*dE3*dN - 4*dN3*dE))
  
  Y <- with (Coeff, 
            D0 + D1*dE + C1*dN + D2*(dE2 - dN2) + C2*(2*dE*dN) +
              D3*(dE3 - 3*dE*dN2) + C3*(3*dE2*dN - dN3) +
              D4*(dE4 - 6*dE2*dN2 + dN4) + C4*(4*dE3*dN - 4*dN3*dE))
  data.frame(X=X, Y=Y)
}
