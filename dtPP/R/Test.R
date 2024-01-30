# Test two models, fit1 and fit2, where fit1 is simple and fit 2 is complex model
# using F-test: fit 2 is significantly better.
# F-test = [(RSS1-RSS2)/(npar2-npar1) ] / [RSS2/(n-np2)]
# http://en.wikipedia.org/wiki/F-test
Ftest.fitPI <- function(fit1, fit2) {

          n <- length(fit1$response)
          np1 <- length(fit1$par)
          np2 <- length(fit2$par)
          Fv <- ( (fit1$ssr - fit2$ssr) / (np2 - np1 ) )/  
                (  fit2$ssr             /  (n - np2 ) )
          
          if (Fv == 0 | is.na(Fv)) 
            p <- 1
          else        
#            p <- 1-pf( Fv, n - np2, n - np1)
            p <- 1-pf( Fv, np2 - np1, n - np2)          

      return(list(F = Fv, p = p, np2 = np2, np1 = np1, n = n))
    }

AIC.fitPI <- function(object, ...) {
  class(object) <- "lm" # quick and dirty
  AIC(object)
}
