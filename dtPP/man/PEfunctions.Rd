\name{PIfunctions}
\alias{PIfunctions}
\alias{fWebb} 
\alias{fWebb.I}
\alias{fJP}
\alias{fJP.I}
\alias{fPG}
\alias{fPG.I}
\alias{fEP}
\alias{fEP.I}

\title{
  Well known photosynthesis-irradiance (P-I) models
}

\description{
  \code{fWebb, fJP, fPG, fEP}, are well-known PI models: Webb (1974); 
  Jassby and Platt (1976); Platt, Galegos and Harrison (1980); Eilers and Peeters (1988);
  
  \code{fWebb.I, fJP.I, fPG.I, fEP.I}, are the standardized versions (respective to I).
}

\usage{
fWebb   (I, p=c(alpha=0.2, ek=200)) 
fWebb.I (I, p=c(alpha=0.2, ek=200)) 

fJP     (I, p=c(alpha=0.15, ek=250))  
fJP.I   (I, p=c(alpha=0.15, ek=250))  

fPG     (I, p=c(alpha=0.2, beta=0, ps=40)) 
fPG.I   (I, p=c(alpha=0.2, beta=0, ps=40))  

fEP     (I, p=c(alpha=0.2, eopt=830, ps=40))  
fEP.I   (I, p=c(alpha=0.2, eopt=830, ps=40))  
}

\arguments{
  \item{p }{The parameter vector or data.frame for the respective P-I model - see details. 
    If a data.frame, then the number of rows should be equal to the length of \code{I}.
  }
  \item{I }{Vector with the photosynthetically active radiation at which the response should be calculated. Default units are microEinst/m2/sec; other units: W/m2. 
  }
}
\value{
 a \code{data.frame} with two columns : the photosynthetically active radiation and the photosynthesis.
}

\details{
 The \code{p} (parameter) arguments for the respective models are (in that order):
 \itemize{
   \item{fWebb or fWebb.I }{: \code{alpha, ek},}
   \item{fJP or fJP.I }{: \code{alpha, ek}} 
   \item{fPG or fPG.I }{: \code{alpha, beta, ps}} 
   \item{fEP or fEP.I }{: \code{alpha, eopt, ps}} 
 } 
 
 where 
 \itemize{
  \item{ \code{ek} and \code{eopt} have the same units as \code{I}; the default: [uEinst/m2/s]}
  \item{ \code{ps} has the same units as \code{pp}, default = [mgC/m3/hr] }
  \item{ \code{alpha}, \code{beta} have units \code{pp/I}, default=[(mgC/m3/hr)/[uEinst/m2/s]}
 }
 
 The models are defined as follows: 
 
 fWebb:  
  \deqn{pp  = alpha*ek*(1-exp(-1*I/ek))}{response  = alpha*ek*(1-exp(-1*I/ek))}

 fWebb.I:
  \deqn{pp  = alpha*ek*(1-exp(-1*I/ek))/I}{response  = alpha*ek*(1-exp(-1*I/ek))/I}

 fJP:
  \deqn{pp = alpha*ek*tanh(I/ek)}{response = alpha*ek*tanh(I/ek)  }

 fJP.I
  \deqn{pp = alpha*ek*tanh(I/ek)/I}{response = alpha*ek*tanh(I/ek)/I}

 fPG:  
  \deqn{pp = ps*(1-exp(-1*alpha*I/ps))*exp(-1*beta*I/ps)}{response = ps*(1-exp(-1*alpha*I/ps))*exp(-1*beta*I/ps) }

 fPG.I:
  \deqn{pp = ps*(1-exp(-1*alpha*I/ps))*exp(-1*beta*I/ps)/I}{response = ps*(1-exp(-1*alpha*I/ps))*exp(-1*beta*I/ps)/I }

 fEP:
  \deqn{pp = I/((1/(alpha*eopt^2))*I^2+(1/ps-2/(alpha*eopt))*I+(1/alpha))}{response = I/((1/(alpha*eopt^2))*I^2+(1/ps-2/(alpha*eopt))*I+(1/alpha))}

 fEP.I
  \deqn{pp = 1/((1/(alpha*eopt^2))*I^2+(1/ps-2/(alpha*eopt))*I+(1/alpha))}{response = 1/((1/(alpha*eopt^2))*I^2+(1/ps-2/(alpha*eopt))*I+(1/alpha))}
 
}

\examples{

## =============================================================
## The different models, default parameter values
## =============================================================

par(mfrow = c(2, 1))
I  <- 0:800

PIcurves <- cbind(
         fWebb(I=I), 
         fJP  (I=I)[,2], 
         fPG  (I=I)[,2], 
         fEP  (I=I)[,2]
         )

matplot(x=PIcurves[,1], y=PIcurves[,-1], type="l", 
        col=1:4, lwd=2, lty=1, 
        xlab="I, uEinst/m2/s", ylab="mgC/m3/h", 
        main="Photosynthesis")
        
legend("bottomright", lty=1, lwd=2, col=1:4, 
       legend=c("Webb", "JP", "PG", "EP"))
       
PIcurves.I <- cbind(
         fWebb.I (I=I), 
         fJP.I   (I=I)[,2], 
         fPG.I   (I=I)[,2], 
         fEP.I   (I=I)[,2]
         )
         
matplot(x=PIcurves.I[,1], y=PIcurves.I[,-1], type="l", 
        col=1:4, lwd=2, lty=1,
        xlab="I, uEinst/m2/s", ylab="mgC/m3/h/[uEinst/m2/s]", 
        main="Normalized photosynthesis")

legend("topright", lty=1, lwd=2, col=1:4, 
       legend=c("Webb.I", "JP.I", "PG.I", "EP.I"))

## =============================================================
## The different models, one parameter value
## =============================================================

par(mfrow = c(2, 2))
I <- 0:150

plot(fWebb(p = c(alpha = 0.5, ek = 180), I = I), 
  main = "Webb", type = "l")
plot(fJP  (p = c(alpha = 0.5, ek = 180), I = I), 
  main = "JP", type = "l")
plot(fPG  (p = c(alpha = 0.5, beta = 1, ps = 180), I = I), 
  main = "PG", type = "l")
plot(fEP  (p = c(alpha = 0.5, eopt = 180, ps = 50), I = I), 
  main = "EP", type = "l")

## =============================================================
## The different models, several parameter values, one I
## =============================================================

par(mfrow = c(2, 2))
I <- 100
alpha <-  seq(0.5, 0.8, length.out=20)
boxplot(fWebb(p = data.frame(alpha=alpha, ek = 180), 
             I = I)$y, 
             main = "Webb")
boxplot(fJP(p = data.frame(alpha=alpha, ek = 180), 
             I = I)$y,  
             main = "JP")
boxplot(fPG  (p = data.frame(alpha=alpha,beta = 1, ps = 180), 
             I = I)$y, 
             main = "PG")
boxplot(fEP  (p = data.frame(alpha=alpha, eopt = 180, ps = 50), 
             I = I)$y, 
             main = "EP")

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
  \code{\link{fitPI}} for fitting functions.

  \code{\link{FRRFdata}}, \code{\link{LABSTAFdata}} for data sets.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  Gregg Silsbe <gregg.silsbe@nioz.nl>
  Jacco Kromkamp <jacco.kromkamp@nioz.nl>
}
\keyword{utilities}
