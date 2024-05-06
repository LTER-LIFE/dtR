
## ==========================================
## ==========================================
## make datasets globally available
## ==========================================
## ==========================================

if (getRversion() >= "2.15.1")  
  utils::globalVariables(c("Wad_depth", "Wad_sediment", 
                           "Wad_Waterheight_HR", "Wad_watertemp_HR", 
                           "Wad_shape", 
                           "Wad_biogeo_RWS"))
