\name{fitPI}
\alias{fitPI}
\alias{fitmultiPI}
\alias{logLik.fitPI}
\alias{predict.fitPI}
\alias{predict.fitmultiPI}
\alias{summary.fitmultiPI}
\alias{print.summary.fitmultiPI}
\alias{r.squared}
\alias{coef.fitmultiPI}
\alias{plot.fitPI}

\title{
  Fits a photosynthesis-irradiance model to one or more sets of P-I data
}

\description{
  \code{fitPI} and \code{fitmultiPI} generate best-fit parameters of a specific PI model;

  \code{predict}, \code{summary}, \code{plot.fitPI} generate predictions or summaries or plot model-data comparisons;
  
  \code{coef.fitmultiPI} creates a table with fitted parameters and environmental variables; 
  
  \code{logLik.fitPI} calculates the log-likelihood of a fit, allowing to estimate the AIC and BIC.

}

\usage{
fitPI (model, I, response, normalized = FALSE, pini = NULL, 
  modargs = NULL, lower = 0, upper = Inf, checkoutlier = FALSE, 
  min.rsq = 0.995, ...)

fitmultiPI(model, id, I, response, data = NULL, normalized = FALSE, pini = NULL, ...)

r.squared(object)

\method{coef}{fitmultiPI}(object, ...)

\method{predict}{fitPI}(object, I, modargs = NULL, normalized = FALSE, ...)

\method{logLik}{fitPI}(object, ...)

\method{predict}{fitmultiPI}(object, I, ...)

\method{summary}{fitmultiPI}(object, ...)

\method{print}{summary.fitmultiPI}(x, digits = max(3, getOption("digits") - 3), ...)

\method{plot}{fitPI}(x, modargs = NULL, normalized = FALSE, ...)
}

\arguments{
  \item{model }{The name of a model or an actual model \code{function} that fits the 
    model to data; names can be one of 'Webb', 'JP', 'PG' or 'EP' for the 
    Webb et al. (1974) model, the Jassby-Platt, Platt-Galegos and the 
    Eilers-Peeters model respectively. See \link{fWebb}, \link{fWebb.I}. 
    When a function, it should have two arguments, \code{p} and \code{I} 
    with \code{p} the parameter values and \code{I} the (photosynthetic) radiation 
    respectively; it should return a data.frame with two columns, named \code{x} and \code{y},
    containing the light (I) and the modeled response. 
  }
  \item{normalized }{When \code{TRUE}, the model is normalized to irradiance;
    this can be useful for fluorescence measures.
  }
  \item{I }{The (photosynthetically active) radiation at which the measurements are made,
    or at which the response should be calculated, a vector with the data, or the 
    name of the column in \code{data}.

   For \code{getPIdata} it is the name of the (independent) light (I) variable to select from the data; 
    if \code{NULL} then it is created based on a variable called 'I'.
  }
  \item{pini }{Inititial guess of parameter values.
  }
  \item{modargs }{Only used if \code{model} is a function: a \code{list} with 
    additional, named arguments passed to the model. See note, examples
  }
  \item{upper, lower }{Upper and lower bounds on the parameter values. 
  
  }
  \item{response }{The actual response value corresponding to \code{I}. 
    Should be a vector of length equal to the length of \code{I}, or the 
    name of the column in \code{data}.
  }
  \item{id }{The identifiers for each PI data set, either a vector of 
  length equal to the length of \code{I}, or the name of the column in \code{data}.
  }
  \item{data }{A PI data set, which has at least a column named ('dbar'). 
    This will be used to create the light ('I') column (if it exists, it will thus be overwritten!), 
    and will generate identifier column called 'id'.
  }
  \item{checkoutlier, min.rsq }{If \code{checkoutlier} is \code{TRUE} and the
  r squared value of the fit is below argument \code{min.rsq} then the model
  will be fit for all cases with one point removed, and the fit which has lowest 
  AIC (akaike information criterium) will be returned. 
  The point removed will be added to the output (item \code{removed}). 
  }
  \item{object }{A PI fit object, as returned by \code{fitPI} or \code{fitmultiPI}.
  }
  \item{x }{A summary of a PI fit, as returned by \code{summary}.
  }
  \item{digits }{The number of significant digits in the printout. 
  }
  \item{... }{Any arguments passed to the fitting function \link{modFit} (package FME), or to the plotting function. 
  }
}
\value{
  Function \code{fitPI} returns a list of class \code{fitPI} that inherits from \link{modFit} from the R-package \code{FME} (and contains items \code{par}, \code{hessian}, \code{residuals}, \code{info}, \code{message}, \code{iterations}, \code{rsstrace}, \code{ssr}, \code{diag}, \code{ms},                 \code{var_ms_unscaled}, \code{var_ms_unweighted}, \code{var_ms}, \code{rank}, \code{df.residual}). 
  
  This means that all methods defined for class \code{modFit} can also be used, such as 
  \code{summary}, \code{deviance}, \code{coef}, \code{residuals}, \code{df.residuals}, 
  \code{plot}, \code{print}. 
  
  In addition, function \code{fitPI} returns the input \code{I}, \code{response}, \code{normalized}, and the (name of the) \code{model}. 
  
  If \code{checkoutlier} is TRUE, it will also return the item \code{removed}, denoting the point considered an outlier and that has been removed. 


  Function \code{fitmultPI} returns a list of class \code{multFRRF}. 
  The only method that works on this class is \code{predict}.
  
  Function \code{r.squared} returns an estimate of the r2 value (not corrected); can be used to rapidly screen good from bad fits.
}

\details{
 The \code{p} (parameter) arguments for the respective models are (in that order):
 \itemize{
   \item{Webb or Webb.I }{: \code{alpha, ek}}
   \item{JP or JP.I }{: \code{alpha, ek}} 
   \item{PG or PG.I }{: \code{alpha, beta, ps}} 
   \item{EP or EP.I }{: \code{alpha, eopt, ps}} 
 }
 See \link{fWebb}, \link{fWebb.I} for the definitions of these models.
  
 Units of the paramters differ based on the data that were fitted:
 \itemize{
  \item{ \code{ek} and \code{eopt} have the same units as \code{I}; the default: [uEinst/m2/s]}
  \item{ \code{ps} has the same units as \code{pp}, default = [mgC/m3/hr] }
  \item{ \code{alpha}, \code{beta} have units \code{pp/I}, default=[(mgC/m3/hr)/[uEinst/m2/s]}
 }

 
}

\note{
 Function \code{plot.fitPI} will predict responses over a sequence of 100 light values
 from the minimum input \code{I} to maximum of input \code{I}.
 
 This means that, if a new model function is defined that has additional arguments
 (i.e. modargs is not \code{NULL}) then 100 values of these arguments, 
 consistent with the 100 \code{I} values should be passed.
}

\examples{

## =============================================================
## Fit existing models to one data set 
## =============================================================

# the 'data set'
I        <- seq(0, 100, length.out = 11)
response <- c(0, 1.6, 1.85, 1.95, 1.99, 2.0, 1.98, 1.95, 1.92, 1.88, 1.82)

# Fit different models
FitWebb <- fitPI(model    = "Webb", 
                 I        = I, 
                 response = response)
FitPG   <- fitPI(model    = "PG",   
                 I        = I, 
                 response = response)
FitJP   <- fitPI(model    = "JP",   
                 I        = I, 
                 response = response)
FitEP   <- fitPI(model    = "EP",   
                 I        = I, 
                 response = response)

# look at some results

summary(FitWebb)
summary(FitPG)
plot(FitWebb)

# use non-default fitting algorithm - not all give reasonable results!

FitW2  <- fitPI(model    = "Webb", 
                I        = I, 
                response = response, 
                method   = "BFGS")
\dontrun{
summary(FitW2)
}

# SSR and r.squared:

c(deviance(FitWebb),  deviance(FitPG),  deviance(FitJP),  deviance(FitEP))

c(r.squared(FitWebb), r.squared(FitPG), r.squared(FitJP), r.squared(FitEP))

# Show fits

par(mfrow = c(1, 1))
plot(I, response, ylim = c(0, 2))

# I for prediction
Ipred <- seq(0, 120, by = 1)

lines(predict(FitWebb, Ipred),                lwd = 2)
lines(predict(FitPG,   Ipred), col = "red",   lwd = 2)
lines(predict(FitJP,   Ipred), col = "blue",  lwd = 2)
lines(predict(FitEP,   Ipred), col = "green", lwd = 2)

legend("bottomright", 
       legend = c("Webb",   "PG",   "JP",    "EP"), 
       col    = c("black", "red", "blue", "green"), 
       lty = 1, lwd = 2)

# use plot function to look at single fits

par(mfrow = c(2,2))

plot(FitWebb, main = "Webb")
plot(FitPG,   main = "Platt-Galegos")
plot(FitJP,   main = "Jassby-Platt")
plot(FitEP,   main = "Eilers-Peeters")
 
## =============================================================
## Fit newly defined model to one data set 
## =============================================================

Monod <- function(I, p) {
  y <- p[1]*I / (I + p[2])
  list(x = I, y = y)
}

I        <- seq(0, 100, length.out = 11)
response <- c(0, 1.6, 1.85, 1.95, 1.99, 2.0, 1.98, 1.95, 1.92, 1.88, 1.82)

FitMonod <- fitPI(model    = Monod, 
                  I        = I, 
                  response = response, 
                  pini     = c(Rm = 2, ks = 25))

# look at some results
summary(FitMonod)
coef(FitMonod)
plot(FitMonod)

par(mfrow = c(1, 1))
plot(I, response, ylim = c(0, 2))

# light for prediction
Ipred <- seq(0, 120, by = 1)

lines(predict(FitMonod, Ipred), lwd = 2)

## =============================================================
## Improving a fit
## =============================================================

# the 'data set' with one outlier
I        <- seq(0, 100, length.out = 11)
response <- c(0, 1.6, 1.85, 1.95, 1.99, 2.0, 1.98, 1.95, 1.92, 1.5, 1.82)

# Fit different models
FitEP   <- fitPI(model    = "EP",   
                 I        = I, 
                 response = response)

plot(FitEP)
AIC(FitEP)

# Remove the outlier manually
FitEP2  <- fitPI(model    = "EP",   
                 I        = I[-10], 
                 response = response[-10])
AIC(FitEP2)                 
lines(predict(FitEP2, I=1:100), col = "red")   

# Automatic removing the outlier
FitEP3  <- fitPI(model    = "EP",   
                 I        = I, 
                 response = response,
                 checkoutlier = TRUE)
                 
lines(predict(FitEP3, I=1:100), col = "blue")  

FitEP3$removed

## =============================================================
## A more complex model depending also on Temperature 
## (Q10 formulation)
## =============================================================

# There are 3 parameters in this model: Max rate (Rm), 
# the half-saturation light intensity (ks) and the temperature coefficient (Q10)

MonodT <- function(I, p, T) {
  
  y <- p[1]*I / (I + p[2]) * p[3]**((T-10)/10)
  list(x = I, y = y)
}

# input data
I        <- seq(0, 100, length.out = 11)
T        <- seq(15, 10, length.out = 11) 
response <- c(0, 1.6, 1.85, 1.95, 1.99, 2.0, 1.98, 1.95, 1.92, 1.88, 1.82)

FitMonodT <- fitPI(model    = MonodT, 
                   I        = I, 
                   response = response, 
                   pini     = c(Rm = 2, ks = 25, Q10 = 2), 
                   modargs  = list(T = T))

# look at some results
summary(FitMonodT)
coef(FitMonodT)

par(mfrow = c(1, 1))

# For plotting the fit, we need input temperature to be of length (100), 
# as the light to predict the response also has 100 values (I).

plot (FitMonodT, 
      modargs = list(T = seq(15, 10, length.out = 100)), 
      main    = "Monod*f(Temp)")

## =============================================================
## Fit models to multiple data sets at once
## =============================================================

# the data generated with tee labSTAF
 set  <- subset(LABSTAFdata, subset = (! is.na(JVPII)))

# The light intensity is called "E" here
# Fit all at once  
 mFitEP <- fitmultiPI(model      = "EP", 
                      data       = set, 
                      id         = "file", 
                      I          = "E", 
                      response   = "JVPII", 
                      normalized = TRUE)
 
 summary(mFitEP)  # look at summary
 coef(mFitEP)     # extract parameters
 r.squared(mFitEP)
 
# plot all fits at once
 plot(mFitEP)

# plot predicted lines in one plot
 par (mfrow = c(1, 1))   
 Par <- 1:1500 
 predEP <- predict (mFitEP, I = Par)
 
 matplot(Par, predEP, 
         type = "l", lwd = 2, col = 1:10, lty = 1)
 
 legend("topright", legend = 1:10, 
        lwd = 2, col = 1:10, lty = 1)
 
 points(LABSTAFdata[,-1], 
        col = (1:10)[LABSTAFdata[,1]], pch = ".", cex = 4)

## =============================================================
## Fit newly defined model to multiple data sets at once
## =============================================================

# newly defined function
 Monod <- function(I, p) {
  y <- p[1]*I / (I + p[2])
  list(x = I, y = y)
 }

 mFitMonod <- fitmultiPI(model      = Monod, 
                         id         = LABSTAFdata$file, 
                         I          = LABSTAFdata$E, 
                         response   = LABSTAFdata$JVPII, 
                         normalized = TRUE, 
                         pini       = c(pm = 1, ks = 1))
                    
 summary(mFitMonod)

# plot predicted lines   
 Par       <- 1:1500 
 predMonod <- predict (mFitMonod, I = Par)
 matplot(Par, predMonod, 
         type = "l", lwd = 2, col = 1:10, lty = 1)


}
\references{
 Jassby A.D. and T. Platt, 1976. Mathematical formulation of the relationship 
 between photosynthesis and light for phytoplankton. Limnol. Oceanogr. 21:540-547.

 Platt T. and A.D. Jasby, 1976. The relationship between photosynthesis and 
 light for natural assemblages of coastal marine phytoplankton. J. Phycol. 12:421-430.

 Platt T., Gallegos C.L. and W.G. Harrison, 1980. Photoinhibition of photosynthesis 
 in natural assemblages of marine phytoplankton. J. Mar. Res., 38, 687-701
 
 Eilers P.H.C. and J.C.H. Peeters, 1988. A model for the relationship between 
 light intensity and the rate of photosynthesis in phytoplankton. Ecol. Modell., 42, 199-215.

 Webb W.L, Newton M. and D. Starr, 1974. Carbon dioxide exchange of Alnus rubra.
 A mathematical model. Oecologia, Berlin 17, 281-291.
}
\seealso{
  \code{\link{fWebb}} for the various PI functions.

  \code{\link{LABSTAFdata}} for a data set.
  
  \code{\link{FRRFdata}} for fitting with normalized and not normalized data.
  
  \code{\link{modFit}} for options of fitting algorithm
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  Gregg Silsbe <gregg.silsbe@nioz.nl>
  Jacco Kromkamp <jacco.kromkamp@nioz.nl>
}
\keyword{utilities}
