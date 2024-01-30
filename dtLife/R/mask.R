# Create a mask if data are in the wadden or not

mask_bat <- function(bat, shape, EPSG=4326){
  MM <- mask_2D(bat$longitude, bat$latitude, shape=shape, EPSG=EPSG)
#  bat$depth[MM==0] <- NA
  attributes(MM)$processing <- c(attributes(MM)$processing,
                                 paste("masked @ ", Sys.time()))
  MM
}
  
mask_2D <- function(longitude, latitude, shape, EPSG=4326){
  xy <-  expand.grid(longitude, latitude)
  names(xy) <- c("longitude","latitude")
  
  xy_sf <- st_as_sf(xy, coords = c("longitude","latitude"))
  xy_sf <- st_set_crs(xy_sf, EPSG)
  
  #if(EPSG!= 4326){
  #str <- "EPSG:4326" 
  #xy_sf = st_transform(xy_sf, str)
  #}
  
  if (missing(shape))
    stop("shapefile 'shape' not given")
  SS    <- st_intersects(xy_sf, shape)
  mask  <- matrix(nrow=length(longitude),
                  ncol=length(latitude),
                  data = unlist(as.numeric(SS)))
#  mask[is.na(mask)] <- 0
  return(list(longitude=longitude, latitude=latitude, 
                     mask = mask))
}

mask_xy <- function(coordinates, shape, EPSG=4326){

  xy_sf <- st_as_sf(coordinates, coords = names(coordinates))
  xy_sf <- st_set_crs(xy_sf, EPSG)
  
  if (missing(shape))
    stop("shapefile 'shape' not given")
  
  SS    <- st_intersects(xy_sf, shape)
  mask  <- unlist(as.numeric(SS))
  return(data.frame(coordinates, mask = mask))
}
