# ==============================================================================
# ==============================================================================
# Masking a bathymetry or other variable with a shape file
# ==============================================================================
# ==============================================================================

mask_shape <- function(coordinates, 
                       longitude, latitude, 
                       shape, # shape file
                       bat, 
                       EPSG = 4326){
  
  if (! missing(bat)) {
    return(mask_bat(bat, shape, EPSG))
    
  } else if (missing(coordinates)){
    
    if (missing(longitude) | missing(latitude))
      stop("either'coordinates' or 'latitude' and 'longitude' should be given")
    
    return(mask_2D(longitude, latitude, shape, EPSG))
  
  } else {
  
    return(mask_xy(coordinates, shape, EPSG))
  }
  
  
}

# ==============================================================================
# Create a mask for a list (bat) with latitude, longitude, 
# ==============================================================================

mask_bat <- function(bat, shape, EPSG = 4326){
  
  MM <- mask_2D(bat$longitude, bat$latitude, 
                shape = shape, 
                EPSG  = EPSG)

  MM$processing <- c(bat$processing,
                     paste("masked @ ", Sys.time()))

  MM
}
  
# ==============================================================================
# Create a mask for latitude, longitude (vectors)  
# ==============================================================================

mask_2D <- function(longitude, latitude, 
                    shape,          # shapefile
                    EPSG = 4326){   # coordinate system
  
  xy        <-  expand.grid(longitude, latitude)
  names(xy) <- c("longitude", "latitude")
  
  xy_sf <- st_as_sf(xy, 
                    coords = c("longitude", "latitude"))
  xy_sf <- st_set_crs(xy_sf, EPSG)
  
  if (missing(shape))
    stop("shapefile 'shape' not given")
  
  SS    <- st_intersects(xy_sf, shape)
  mask  <- matrix(nrow = length(longitude),
                  ncol = length(latitude),
                  data = unlist(as.numeric(SS)))
#  mask[is.na(mask)] <- 0
  return(list(longitude = longitude, 
              latitude  = latitude, 
              mask      = mask,
              asp       = aspect_coord(latitude)
  ))
}

# ==============================================================================
# Create a mask for a data.frame with coordinates 
# ==============================================================================

mask_xy <- function(coordinates, 
                    shape, 
                    EPSG = 4326){
  
  if (! is.data.frame(coordinates))
    coordinates <- as.data.frame(coordinates)
  
  xy_sf <- st_as_sf(coordinates, 
                    coords = names(coordinates))
  xy_sf <- st_set_crs(xy_sf, EPSG)
  
  if (missing(shape))
    stop("shapefile 'shape' not given")
  
  SS    <- st_intersects(xy_sf, shape)
  mask  <- unlist(as.numeric(SS))
  
  return(data.frame(coordinates, mask = mask))
}
