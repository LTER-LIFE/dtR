
## ==========================================
## ==========================================
## make datasets globally available
## ==========================================
## ==========================================

if (getRversion() >= "2.15.1")  
  utils::globalVariables(c("WadDepth", "WadGrain", 
                           "WadHeightHR", "WadTempHR", 
                           "WadShape", 
                           "WadBiogeo"))
