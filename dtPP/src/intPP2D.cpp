#include <Rcpp.h>
using namespace Rcpp;

// 
// C++ functions to estimate depth-integrated primry production
//

// [[Rcpp::plugins("cpp11")]]

// =============================================================================
// =============================================================================
//
// Integrated primary production (PP), assuming a mixed water column
// Uses 5-points gaussian integration to integrate PP over depth.
//
// The photosynthesis-light relationship is represented with the 
// Eilers-Peeters model (parameters alpha, eopt and pmax).
//
// The extinction coefficients, light intensity, water surface elevation,
// and the photosynthesis parameters change over time
//
// 2-D versions
//
// =============================================================================
// =============================================================================

// [[Rcpp::export]]

NumericMatrix intPP_mixed(      // spatial-temporal integrated production function

       NumericVector   zmax,    // water depth                         (ns)
       NumericMatrix     kd,    // extinction coefficient              (nt, ns)
       NumericMatrix    par,    // photosynthetically active radiation (nt, ns)
       NumericMatrix   alfa,    // slope of PI curve at par=0          (nt, ns)
       NumericMatrix   eopt,    // optimal light intensity             (nt, ns)
       NumericMatrix   pmax,    // maximum productivity                (nt, ns)
       NumericMatrix height) {  // water surface elevation             (nt, ns)

// ========================
// declaration section
// ========================
  
  double parz , 
         zt, zg, zz,
         total,
         f1, f2, f3,
         ps   ;
  
  // abscissas of 5 gaussian points 
  static double xg[] = {0.04691008, 0.23076534,    0.5,   0.76923466, 0.95308992};
  
  // weight of the points 
  static double wg[] = {0.1184634,  0.2393143, 0.2844444, 0.2393143,  0.1184634};
  
  int  j, it, is;            // iterators
  int  nt = kd.nrow();       // sizes
  int  ns = zmax.size();
  
  // create a matrix filled with 0's (for output)
  NumericMatrix pp(nt, ns);
  
// ========================
// estimate photosynthesis 
// ========================
  
  // loop over all stations
  
  for (is = 0; is < ns; ++is) {
 
  // loop over time
    for (it = 0; it < nt; ++it) {              
    
      zt  = zmax[is] + height(it, is);      // current water depth
      
        
      if (zt > 0.0 && (par(it, is) > 0)) {           // there is water -> PP
      
      // depth where light = 0.1 % of total (euphotic depth) or where bottom
        zg    = -log(0.001)/kd(it, is); 
        if (zg > zt) zg = zt;
        
        // elements of the Eilers-Peeters model for this station
        f1  = 1.0/(alfa(it, is) * pow(eopt(it, is), 2.0));
        f2  = 1.0/pmax(it, is) - 2.0/(alfa(it, is) * eopt(it, is));
        f3  = 1.0/alfa(it, is);
        
        total = 0.0;
    
      // loop over 5 gaussian points
        for (j = 0; j < 5; ++j) {              
          zz   = xg[j]*zg;                 // depth at which to estimate ps
        
          parz = par(it, is) * exp (-kd(it, is) * zz); // light at this depth 
        
          // photosynthesis by the eilers-peeters model
          ps = parz/(f1 * pow(parz, 2.0) + f2 * parz + f3);          
          
          total += wg[j]*ps*zg;            // integrated pp
        }
      } else total = 0.0;
    
      pp(it, is) = total;                   // save integrated pp
    }
  }
  // ========================
  // results
  // ========================
  
  return pp;
}


// =============================================================================
// =============================================================================
//
// Integrated primary production (PP), assuming a mixed water column
// Uses 5-points gaussian integration to integrate PP over depth.
//
// The photosynthesis-light relationship is represented with the 
// Eilers-Peeters model (parameters alpha, eopt and pmax).
//
// The extinction coefficients, light intensity, water surface elevation,
// and the photosynthesis parameters all change over time
//
// =============================================================================
// =============================================================================

// [[Rcpp::export]]

NumericMatrix rad_bot(       // radiation at the bottom of the water column
    
    NumericVector   zmax,    // water depth                         (ns)
    NumericMatrix     kd,    // extinction coefficient              (nt, ns)
    NumericMatrix    par,    // photosynthetically active radiation (nt, ns)
    NumericMatrix height) {  // water surface elevation             (nt, ns)
  
  // ========================
  // declaration section
  // ========================
  
  int  it, is;            // iterators
  int  nt = kd.nrow();
  int  ns = zmax.size();
  
  double zt; 
  
  // create a vector filled with 0's (for output)
  NumericMatrix rad(nt, ns);
  
  // ========================
  // estimate photosynthesis 
  // ========================
  
  // loop over all stations
  for (is = 0; is < ns; ++is) {
    
    // loop over time
    for (it = 0; it < nt; ++it) {              
      
      zt  = zmax[is] + height(it,is);      // current water depth
      
      if (zt > 0.0) {                      // water on top of sediment
              rad(it,is) = par(it,is) * exp (-kd(it,is) * zt);
      } else  rad(it,is) = par(it,is)  ;  
    }
  }
  
  // ========================
  // results
  // ========================
  
  return rad;
}

// =============================================================================
// =============================================================================
//
// Integrated primary production (PP), assuming an exponentially decreasing
// chlorophyll concentration.
//
// Uses 5-points gaussian integration to integrate PP over depth.
//
// The photosynthesis-light relationship is represented with the 
// Eilers-Peeters model (parameters alpha, eopt and pmax), 
// where alpha and pmax are chlorophyll-specific values.
// 
// For instance, pmax has units of mgC/mgChl/h. 
//
// The chlorophyll concentration is estimated based on the mud content
//
// The extinction coefficients, light intensity, 
// and the photosynthesis parameters all change over time
//
// =============================================================================
// =============================================================================


// [[Rcpp::export]]
NumericMatrix intPP_exp(       
        NumericVector     zmax,  // sediment depth                        (ns)
        NumericVector       kd,  // extinction coefficient                (ns)
        NumericVector     pMud,  // fraction of mud in sediment           (ns)        
        NumericMatrix      par,  // photosynthetically active radiation (nt, ns)
        NumericMatrix     alfa,  // slope of PI curve                   (nt, ns)
        NumericMatrix     eopt,  // optimal light intensity             (nt, ns)
        NumericMatrix     pmax){ // maximum productivity                (nt, ns)
                                
  // ========================
  // declaration section
  // ========================
  
  double parz , zt, zg, zz,
         total, alfa_z, pmax_z, eopt_z, pChl_z, 
         f1, f2, f3,
         ps   ;
  
  // abscissas of 5 gaussian points 
  static double xg[] = {0.04691008, 0.23076534,    0.5,   0.76923466, 0.95308992};
  
  // weight of the points 
  static double wg[] = {0.1184634,  0.2393143, 0.2844444, 0.2393143,  0.1184634};
  
  int  j, it, is;            // iterators
  int  nt = par.nrow();
  int  ns = zmax.size();
  
  // create a vector filled with 0's (for output)
  NumericMatrix pp(nt, ns);
  
  // ========================
  // estimate photosynthesis 
  // ========================
  
  // loop over all stations
  
  for (is = 0; is < ns; ++is) {
    
    
    // loop over time
    for (it = 0; it < nt; ++it) {              
      
      zt  = zmax[is] ;          // sediment depth where there is chlorophyll

      if (zt > 0.0 && (par(it,is) > 0)) {           // there is sediment -> PP
        
        
        // depth where light = 0.1 % of total (euphotic depth) or where no more Chl
        zg    =  -log(0.001)/kd[is];
        if (zg > zt) zg = zt; 
        
        // loop over 5 gaussian points (depth)
        total = 0.0;
        
        for (j = 0; j < 5; ++j) {      
          
          // depth at which to estimate ps
          
          zz    = xg[j]*zg;                 
          
          // Chl/Chlmean at this depth
          pChl_z  = (1-pMud[is]) + pMud[is]*2.0*exp(1.0)*exp(-2.*exp(1.0)*zz/zt);
          
          // Eilers-Peeters parameters for this chlorophyll conc
          alfa_z = alfa(it, is) * pChl_z;  
          pmax_z = pmax(it, is) * pChl_z;
          eopt_z = eopt(it, is);
          
          // elements of the Eilers-Peeters model for this station
          f1  = 1.0/(alfa_z * pow(eopt_z, 2.0));
          f2  = 1.0/pmax_z - 2.0/(alfa_z * eopt_z);
          f3  = 1.0/alfa_z;
          
          // light at this depth  
          parz = par(it,is) * exp (-kd[is] * zz); 
          
          // photosynthesis by the eilers-peeters model
          ps = parz/(f1 * pow(parz, 2.0) + f2 * parz + f3);          
          
          total += wg[j]*ps*zg;            // integrated pp
        }
      } else total = 0.0;
      
      pp(it,is) = total;                   // save integrated pp
    }
  }
  // ========================
  // results
  // ========================
  
  return pp;
}



// =============================================================================
// =============================================================================
//
// Integrated primary production (PP), assuming a normally distributed
// chlorophyll concentration.
//
// Uses 5-points gaussian integration to estimated integrated chlorophyll
// and to integrate PP over depth.
//
// The photosynthesis-light relationship is represented with the 
// Eilers-Peeters model (parameters alpha, eopt and pmax), 
// where alpha and pmax are chlorophyll-specific values.
// 
// For instance, pmax has units of mgC/mgChl/h. 
//
// The chlorophyll concentration is estimated based on the mud content
//
// The extinction coefficients, light intensity, 
// and the photosynthesis parameters all change over time
//
// =============================================================================
// =============================================================================


// [[Rcpp::export]]
NumericMatrix intPP_log(       
    NumericVector     zmax,  // sediment depth                        (ns)
    NumericVector       kd,  // extinction coefficient                (ns)
    NumericVector       zn,  // depth where chl peaks                 (ns)        
    NumericVector      zsd,  // sd of chl peak                        (ns)        
    NumericMatrix      par,  // photosynthetically active radiation (nt, ns)
    NumericMatrix     alfa,  // slope of PI curve                   (nt, ns)
    NumericMatrix     eopt,  // optimal light intensity             (nt, ns)
    NumericMatrix     pmax){ // maximum productivity                (nt, ns)
  
  // ========================
  // declaration section
  // ========================
  
  double parz , zt, zg, zz,
  total, alfa_z, pmax_z, eopt_z, fac, pChl_z,
  f1, f2, f3,
  ps   ;
  
  // abscissas of 5 gaussian points 
  static double xg[] = {0.04691008, 0.23076534,    0.5,   0.76923466, 0.95308992};
  
  // weight of the points 
  static double wg[] = {0.1184634,  0.2393143, 0.2844444, 0.2393143,  0.1184634};
  
  int  j, it, is;            // iterators
  int  nt = par.nrow();
  int  ns = zmax.size();
  
  // create a vector filled with 0's (for output)
  NumericMatrix pp(nt, ns);
  
  // ========================
  // estimate photosynthesis 
  // ========================
  
  // loop over all stations
  for (is = 0; is < ns; ++is) {
    
    // loop over time
    for (it = 0; it < nt; ++it) {              
      
      // depth above which 99.9% of chlorophyll
      zt  = R::qnorm(0.999, zn[is], zsd[is], true, false);
      
      if (zt > zmax[is])
        zt  = zmax[is] ;          // depth over which to integrate
      
      if (zt > 0.0 && (par(it,is) > 0)) {  // there is chlorophyll AND light

        total = 0.0;
        
        // loop over 5 gaussian points (depth) to estimate chlorophyll
        for (j = 0; j < 5; ++j) {      
          zz    = xg[j]*zt;      // zt = depth where there is chlorophyll                 
          
          // relative Chl at this depth: normal distr pts at zz, mean zn, sd zsd
          pChl_z  =  R::dnorm(zz, zn[is], zsd[is], false);
          
          total  += wg[j] * pChl_z * zt;
        }
        
        fac = 1.0/total;
        
        // loop over 5 gaussian points (depth) to estimate photosynthesis
        total = 0.0;
        
        // depth where light = 0.1 % of total (euphotic depth) or where no more Chl
        zg    =  -log(0.001)/kd[is];
        if (zg > zt) zg = zt; 
        
        for (j = 0; j < 5; ++j) {      
          
          zz    = xg[j]*zg;        // zg = depth where there is stil 1% light         
          
          // Chl/Chlmean at this depth
          pChl_z  = fac * R::dnorm(zz, zn[is], zsd[is], false);
          
          // Eilers-Peeters parameters for this chlorophyll conc
          alfa_z = alfa(it, is) * pChl_z;  
          pmax_z = pmax(it, is) * pChl_z;
          eopt_z = eopt(it, is);
          
          // elements of the Eilers-Peeters model for this station
          f1  = 1.0/(alfa_z * pow(eopt_z, 2.0));
          f2  = 1.0/pmax_z - 2.0/(alfa_z * eopt_z);
          f3  = 1.0/alfa_z;
          
          // light at this depth  
          parz = par(it, is) * exp (-kd[is] * zz); 
          
          // photosynthesis by the eilers-peeters model
          ps = parz/(f1 * pow(parz, 2.0) + f2 * parz + f3);          
          
          total += wg[j]*ps*zg;            // integrated pp
        }
      } else total = 0.0;
      
      pp(it, is) = total;                   // save integrated pp
    }
  }
  return(pp);
}
