test_that("simple interpolate_xyt works", {

  times <- seq(1, 10, length.out=1000)

  d1   <- cbind(x     = 1, 
                y     = 2,   
                times = times, 
                value = sin(times))
  d2   <- cbind(x     = 1, 
                y     = 3,   
                times = times, 
                value = sin(times+0.5))
  d3   <- cbind(x     = 2, 
                y     = 2,   
                times = times, 
                value = sin(times-0.5))
  d4   <- cbind(x     = 2, 
                y     = 3,   
                times = times, 
                value = sin(times-1))
  d5   <- cbind(x     = 2.5, 
                y     = 2.5, 
                times = times, 
                value = sin(times+1))
  
  # combine all to one input data set
  input_xytv <- rbind(d1, d2, d3, d4, d5)
  
  # 4 points and 10 time values to interpolate to
  
  output_xy <- cbind(x = runif(n=4, min=1, max=2), 
                     y = runif(n=4, min=2, max=3))
  output_t  <- seq(4, 8, length.out=10)
  
  # Interpolate to 2D grid
    expect_no_error(output_v <- interpolate_xyt(input_xytv = input_xytv, 
                              output_xy  = output_xy, 
                              output_t   = output_t  ))
})


test_that("2D interpolation of weather data - inverse distance", {
  Weather.stations <- attributes(Wad_weather)$stations

  # Create daily averages of the Windspeed as input
  
  Wad_weather$Day <- as.Date(Wad_weather$datetime)
  
  mf <- par(mfrow=c(3,1))
  
  WS <- aggregate(x   = Wad_weather$windspeed, 
                  by  = list(Wad_weather$Day, Wad_weather$station), 
                  FUN = mean)
  colnames(WS) <- c("date", "station", "windspeed")               
  
  Input <- merge(
           Weather.stations[ , c("station", "longitude", "latitude")],
           WS              [ , c("station", "date", "windspeed")]
                )
                
  # remove first column                 
  Input     <- Input[, -1] 
  
  # head(Input)   # 4 columns: x, y, t, value
  
  # output: 400 stations, first day of the month
  
  # x-and y values to output
  nx   <- 50
  ny   <- 50
  
  outx <- seq(4.8,   5.5, length.out=nx)
  outy <- seq(52.9, 53.4, length.out=ny)
  
  output_xy <- expand.grid(x= outx,
                           y= outy)
  head(output_xy)
  
  # output time value
  output_t <- as.Date(c("2021-01-01", "2021-02-01", 
                        "2021-03-01", "2021-04-01", 
                        "2021-05-01", "2021-06-01", 
                        "2021-07-01", "2021-08-01", 
                        "2021-09-01", "2021-10-01", 
                        "2021-11-01", "2021-12-01"))

  expect_no_error(output_wind <- interpolate_xyt(input_xytv  = Input, 
                               output_xy   = output_xy, 
                               barycentric = FALSE,
                               output_t    = output_t))
  
})