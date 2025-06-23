// =============================================================================
// =============================================================================
//
// Interpolation in space using different input and output formats
// Uses weighted interpolation, with inverse distance weighing
// author: Karline soetaert
//
// =============================================================================
// =============================================================================


#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

// =============================================================================
// =============================================================================
//
// accessory functions
//
// =============================================================================
// =============================================================================

// =============================================================================
// find maximum value in a small vector
// =============================================================================

int update_imax(double *w, 
                int nmean){
  
  int imax = 0;
  
  for (int iw = 1; iw < nmean; iw++){
    
    if (w[iw] > w[imax]) imax = iw;
      
  }
  return imax;
}

// =============================================================================
// update weight and value with smaller value (if smaller)
// =============================================================================

void update_w_v(double *w, double ww, 
                double *v, double vv, 
                int *imax,
                int ii, int nmean) {
  
  if (ii < nmean) {  // vectors w and v not yet filled
    
    w[ii] = ww;
    v[ii] = vv;
    
  } else {              
    
    if (ii == nmean)      // find  maximal value (first to replace)
      *imax = update_imax(w, nmean);
    
    if (ww < w[*imax]) {  // check if smaller
      
      w[*imax] = ww;      // replace
      v[*imax] = vv;
      *imax = update_imax(w, nmean);
    }  
    
  } // end else
}

// ============================================================================
// estimate weights from squared distance 
// - check for exact matches
// ============================================================================

void inv_distance (double *w,  // input: distance^2, output: weight
                   int nmean){ // number of elements in w     
  
  bool   hasNAN = false;
  int    iw;
  double sum_w;
  
  for (iw = 0; iw < nmean; ++iw) { 
    
    if (w[iw] == 0) {          // exact match 
      
      hasNAN = true;
      w[iw] = -99.;            // negative values are easy to find
      
    } else  
      w [iw] = 1./sqrt(w[iw]); // from distance^2 to weight
  } 
  
  if (hasNAN){                 // there are exact matches
    
    for (iw = 0; iw < nmean; ++iw) { 
      
      if (w[iw] < 0 ) {        // exact match
        w[iw] = 1;
      } else w[iw] = 0;
    }
    
  }  
  
  // rescale weights    
  sum_w = 0.;
  
  for (iw = 0; iw < nmean; ++iw) { 
    sum_w  += w[iw];
  } 
  
  if (sum_w > 0){
    
    for (iw = 0; iw < nmean; ++iw) { 
      w[iw] = w[iw] / sum_w;
    }
    
  }
}

// ============================================================================
// interpolate 
// ============================================================================

double get_value (double *w,  // input: distance^2, output: weight
                  double *v,  // values - can be NA
                  int nmean){ // number of elements in w     
  
  bool   hasNAN = false;
  int    iw;
  double sum_w;
  double res;
  
  // estimate weights

  for (iw = 0; iw < nmean; ++iw) { 
    
    if (w[iw] == 0) {  // exact match 
      
      hasNAN = true;
      w[iw]  = -99.;   // negative values are simple to find
    } else if (w[iw] > 0.9e100) { // point was too far
      w[iw] = 0.;  
    } else  
      w [iw] = 1./sqrt(w[iw]); // from distance^2 to weight
  } 
  
  // check for NA in inputs
    
  for (iw = 0; iw < nmean; ++iw) { 
    
    if (v[iw] == NA_REAL){
      w[iw] = 0;    // set weight = 0
      v[iw] = 0;    // value not important, but not NA
    
    } 
  } 
  
  if (hasNAN) {     // there are exact matches
    
    for (iw = 0; iw < nmean; ++iw) { 
      
      if (w[iw] < 0 ) {
        w[iw] = 1;      // exact match gets weight = 1
      } else w[iw] = 0; // no    match gets weight = 0
    }
    
  }  
  
  // rescale weights    
  sum_w = 0.;
  
  for (iw = 0; iw < nmean; ++iw) { 
    sum_w  += w[iw];
  } 
  
  // weighted average
  
  
  if (sum_w > 0) {
    // rescale w
    for (iw = 0; iw < nmean; ++iw) { 
      w[iw] = w[iw] / sum_w;
    }
    
    // weighted average
    res = 0.;
    for (iw = 0; iw < nmean; ++iw) { 
      res += w[iw]*v[iw];
    }
    
  } else res = NA_REAL;  // sum_w = 0 if all values are NAs
  
  return res;
}

// =============================================================================
// =============================================================================
//
// input : 2D data, 
// output: xy format
//
// interpolated over 4 corners in the grid cell where the output point falls
//
// =============================================================================
// =============================================================================

// [[Rcpp::export]]

NumericVector interpolate_2D_xy_cpp(  
    
    NumericVector   input_x,  //  (nix)
    NumericVector   input_y,  //  (niy)
    NumericMatrix  input_2D,  //  (nix, niy)
    NumericMatrix output_xy,  //  (nout, 2)
    double             asp )  //  y/x aspect ratio  
  
{  
  
  // ------------------------
  // declaration section
  // ------------------------
  
  double x, y, dd            ;
  double w[4], v[4]   ;           // weights and values
  
  int  nix =    input_x.size();   // input/output sizes
  int  niy =    input_y.size();       
  int  nout = output_xy.nrow();       
  
  int  ix, iy, io;                // iterators
  
  NumericVector vv(nout);         // output
  
  // ------------------------
  // calculation 
  // ------------------------
  
  // loop over all required outputs
  // will use 4 points to average, stored as:
  // 1 3
  // 0 2
  
  for (io = 0; io < nout; ++io) {    // output points
    
    x = output_xy(io, 0);
    if (x > input_x[nix-1]) x = input_x[nix-1];
    if (x < input_x[0])     x = input_x[0];
    
    y = output_xy(io, 1);
    if (y > input_y[niy-1]) y = input_y[niy-1];
    if (y < input_y[0])     y = input_y[0];
    
    // find grid position in x
    for (ix = 1; ix < nix; ix++) {  // loop over input_x (should be sorted)
      if (input_x[ix] >= x) break; 
    }  
    
    // squared distance in x 
    w[0] = pow(x - input_x[ix-1], 2); // w0 : (x0, y0)
    w[1] = w[0];                      // w1 : (x0, y1)
    w[2] = pow(x - input_x[ix  ], 2); // w2 : (x1, y0)
    w[3] = w[2];                      // w2 : (x1, y1)  
    
    // find grid position iny
    for (iy = 1; iy < niy; iy++) {              
      if (input_y[iy] >= y) break; 
    }  
    
    // add squared distance in y
    dd    = pow(y - input_y[iy-1], 2)*asp*asp; // account for aspect ratio 
    w[0] += dd; 
    w[2] += dd;
    
    dd    = pow(y - input_y[iy  ], 2)*asp*asp;
    w[1] += dd;
    w[3] += dd;
    
    // values at the corners to interpolate
    v[0] = input_2D(ix-1, iy-1) ;
    v[1] = input_2D(ix-1, iy  ) ;
    v[2] = input_2D(ix  , iy-1) ;
    v[3] = input_2D(ix  , iy  ) ;
    
    // weighted averages
    vv(io) = get_value (w, v, 4);
    
  }
  
  return vv;
}

// =============================================================================
// =============================================================================
//
// input : 2D data, 
// output: xy format - regular x and y input
//
// =============================================================================
// =============================================================================

// [[Rcpp::export]]

NumericVector interpolate_2D_xy_r_cpp(  
    
    NumericVector   range_x,  //             (3)
    NumericVector   range_y,  //             (3)
    NumericMatrix  input_2D,  //             (nix, niy)
    NumericMatrix output_xy,  //             (nout, 2)
    double             asp )  // y/x aspect ratio  
  
{  
  
  // ------------------------
  // declaration section
  // ------------------------
  
  double xo, yo, xi, yi, dd  ;
  double w[4]     ;   // weights
  double v[4]     ;   // values
  
  int  io;                   // iterators
  
  // output grid in x
  double xmin = range_x(0);  
  double dx   = range_x(2);
  
  // output grid in y
  double ymin = range_y(0);
  double dy   = range_y(2);
  
  int  nout  = output_xy.nrow();       
  int  nix   = input_2D.nrow();       
  int  niy   = input_2D.ncol();       
  
  NumericVector vv(nout);       // output vector
  
  // ------------------------
  // calculation 
  // ------------------------
  
  for (io = 0; io < nout; ++io) {  // loop over all outputs

    xo = output_xy(io, 0);
    
    // position in input grid, x-direction
    int ix = (xo - xmin)/dx;
    
    if (ix <= 0)   ix = 1;
    if (ix >= nix) ix = nix-1;
    
    xi = (xmin + double(ix-1) * dx ); // min x value to use
    
    yo = output_xy(io, 1);
    
    int iy = (yo - ymin)/dy;
    if (iy <= 0)   iy = 1;
    if (iy >= niy) iy = niy-1;
    
    yi = (ymin + double(iy-1)*dy );
    
    // squared distance in x-direction
    w[0] = pow(xo - xi, 2); 
    w[1] = w[0]; 
    w[2] = pow(xo - (xi + dx), 2);
    w[3] = w[2];
    
    // add squared distance in y-take into account aspect ratio
    dd   = pow(yo - yi, 2)*asp*asp;
    w[0] += dd; 
    w[2] += dd; 
    dd = pow(yo - (yi + dy), 2)*asp*asp;
    w[1] += dd;
    w[3] += dd;
    
    // values at the corners to interpolate
    v[0] = input_2D(ix-1, iy-1) ;
    v[1] = input_2D(ix-1, iy  ) ;
    v[2] = input_2D(ix  , iy-1) ;
    v[3] = input_2D(ix  , iy  ) ;

    // weighted averages
    vv[io] = get_value (w, v, 4);
        
  }
  
  return vv;
}


// =============================================================================
// =============================================================================
//
// input: xy data (long format)
// output: 2D format (wide)
//
// =============================================================================
// =============================================================================

// [[Rcpp::export]]

NumericMatrix interpolate_xy_2D_cpp(  
    
    NumericMatrix input_xyv,  //             (nin, 3)
    NumericVector  output_x,  //             (nox)
    NumericVector  output_y,  //             (noy)
    int               nmean,  // number of points to average 
    double      max_distance,  // above this value: assume NA, unless negative
    double             asp )  // y/x aspect ratio  
  
{  
  
  // ------------------------
  // declaration section
  // ------------------------
  
  double x, y, dd            ;
  double *w = new double[nmean]   ;   // weights
  double *v = new double[nmean]   ;   // values to weigh
  
  int  iox, ioy, iin;            // iterators
  
  int  nin =  input_xyv.nrow();  // sizes
  int  nox =   output_x.size();       
  int  noy =   output_y.size();       
  
  int imax = 0;
  
  double  *D_x = new double[nin];
  double  D_y;
  
  NumericMatrix  vv(nox, noy);      // output
  
  // ------------------------
  // calculation 
  // ------------------------
  
    for (iox = 0; iox < nox; ++iox) {  // output x
      
      // squared distance in x-direction for all iin
      x = output_x(iox);
      
      for (iin = 0; iin < nin; ++iin) {  // input
        D_x[iin] = pow(input_xyv(iin, 0) - x, 2); 
      }  

      for (ioy = 0; ioy < noy; ++ioy) {  // output y
      
        for (iin = 0; iin < nin; ++iin) {  // input
          
          y   = input_xyv(iin, 1);
          
          // squared distance in y-direction 
          D_y = pow(y - output_y(ioy), 2)*asp*asp; 
          
          // squared total distance
          dd = D_x[iin] + D_y;  
          
          if (max_distance > 0){ 
             if (dd > max_distance) dd =  1.0e100;
          }
          
          // update weights w, with dd; value v with input(,2)
          update_w_v( w, dd,                
                      v, input_xyv(iin, 2), 
                      &imax, iin, nmean);
         } // iin  
      
      // weighted averages
      vv(iox, ioy) = get_value (w, v, nmean);
      
    }  // ioy
  }    // iox
  
  return vv;
}


// =============================================================================
// =============================================================================
//
// input: 2D data, 
// output: 2D format
//
// =============================================================================
// =============================================================================

// [[Rcpp::export]]

NumericMatrix interpolate_2D_2D_cpp(  
    
    NumericVector   input_x,  //             (nix)
    NumericVector   input_y,  //             (niy)
    NumericMatrix  input_2D,  //             (nix, niy)
    NumericVector  output_x,  //             (nox)
    NumericVector  output_y,  //             (noy)
    double             asp )  // y/x aspect ratio  
  
{  
  
  // ------------------------
  // declaration section
  // ------------------------
  
  double dd                 ;
  double w[4], v[4]  ;        // weights and values
  
  int  iix, iiy, iox, ioy, iw, ii;   // iterators
  int imax = 0;
  
  int  nix =   input_x.size();       // sizes
  int  niy =   input_y.size();       
  int  nox =  output_x.size();       
  int  noy =  output_y.size();       

  // Set element x1,y1 to 5
  
  double  *D_x = new double[nix];
  double  D_y;
  
  NumericMatrix   vv(nox, noy);      // output
  
  // ------------------------
  // calculation 
  // ------------------------
  
  int nmean = 4;
  
  // loop over output
  for (iox = 0; iox < nox; ++iox) {
    
    for (iix = 0; iix < nix; ++iix) {
        D_x[iix] =  pow(input_x[iix] - output_x(iox), 2); 
    }  
    
    for (ioy = 0; ioy < noy; ++ioy) {
      
      for (iw = 0; iw < nmean; ++iw) {
        w[iw] = 0;
        v[iw] = 0;
      }
      
      ii = 0;
      
      for (iix = 0; iix < nix; ++iix) {
        
        for (iiy = 0; iiy < niy; ++iiy) {
          D_y = pow(input_y[iiy] - output_y(ioy), 2)*asp*asp; 
          dd = D_x[iix] + D_y;  // squared total distance
          
          update_w_v( w, dd, 
                      v, input_2D(iix, iiy), 
                      &imax, ii, nmean);
                      ii++;
                      
        } // input y
      } // input x
      
      // weighted averages
      vv(iox, ioy) = get_value (w, v, nmean);
      
    }  // ioy
  }    // iox
  
  return vv;
}

// same but regular inputs
// [[Rcpp::export]]

NumericMatrix interpolate_2D_2D_r_cpp(  
    
    NumericVector   range_x,  //             (3)
    NumericVector   range_y,  //             (3)
    NumericMatrix  input_2D,  //             (nix, niy)
    NumericVector  output_x,  //             (noy)
    NumericVector  output_y,  //             (noy)
    double             asp )  // y/x aspect ratio  
  
{  
  
  // ------------------------
  // declaration section
  // ------------------------
  
  double xo, yo, xi, yi, dd  ;
  double w[4], wx[4]   ;   // weights
  double v[4]          ;   // values
  
  int  iox, ioy;                   // iterators
  
  // output grid in x
  double xmin = range_x(0);  
  double dx   = range_x(2);
  
  // output grid in y
  double ymin = range_y(0);
  double dy   = range_y(2);
  
  int  nix   = input_2D.nrow();       
  int  niy   = input_2D.ncol();       
  
  int  nox =  output_x.size();       
  int  noy =  output_y.size();       
  
  NumericMatrix   vv(nox, noy);      // output
  
  // ------------------------
  // calculation 
  // ------------------------
  
  for (iox = 0; iox < nox; ++iox) {
    
    xo = output_x(iox);
    
    // position in input grid, x-direction
    int ix = (xo - xmin)/dx;
    if (ix <= 0)   ix= 1;
    if (ix >= nix) ix = nix-1;
    
    xi = (xmin + double(ix-1) * dx ); // min x value to use
    
    // squared distance in x-direction
    wx[0] = pow(xo - xi, 2); 
    wx[1] = wx[0]; 
    wx[2] = pow(xo - (xi + dx), 2);
    wx[3] = wx[2];
      
    for (ioy = 0; ioy < noy; ++ioy) {
      
      yo = output_y(ioy);
      
      int iy = (yo - ymin)/dy;
      if (iy <= 0)   iy = 1;
      if (iy >= niy) iy = niy-1;
      
      yi = (ymin + double(iy-1)*dy );
      
      // add squared distance in y-take into account aspect ratio
      dd   = pow(yo - yi, 2)*asp*asp;
      w[0] = wx[0] + dd; 
      w[2] = wx[2] + dd; 
      dd = pow(yo - (yi + dy), 2)*asp*asp;
      w[1] = wx[1] + dd;
      w[3] = wx[3] + dd;
      
      // values at the corners to interpolate
      v[0] = input_2D(ix-1, iy-1) ;
      v[1] = input_2D(ix-1, iy  ) ;
      v[2] = input_2D(ix  , iy-1) ;
      v[3] = input_2D(ix  , iy  ) ;
      
      
      // weighted averages
      vv(iox, ioy) = get_value (w, v, 4);
      
    }
  }
  return vv;
}
  

// =============================================================================
// =============================================================================
//
// input and output in xy format
//
// =============================================================================
// =============================================================================

// [[Rcpp::export]]

NumericVector interpolate_xy_xy_cpp(  
    
    NumericMatrix input_xyv,  //  (nin, 3)
    NumericMatrix output_xy,  //  (nout, 2)
    int               nmean,  //  number to average  
    double      max_distance,  // above this value: assume NA, unless negative
    double             asp )  // y/x aspect ratio  
  
{  
  // ------------------------
  // declaration section
  // ------------------------
  
  double x, y, dd      ; 
  
  double *w = new double[nmean]   ;   // weights
  double *v = new double[nmean]   ;   // values to weigh
//  double w[nmean], v[nmean] ; // weights, values
  
  int  nin = input_xyv.nrow();       // sizes
  int  nout = output_xy.nrow();       
  
  int  ii, io;                       // iterators
  int  imax = 0;
  
  NumericVector vv(nout);            // output vector
  
  // ------------------------
  // calculation 
  // ------------------------
  
  for (io = 0; io < nout; ++io) {    // loop over output
    
    x = output_xy(io, 0);
    y = output_xy(io, 1);
    
    for (ii = 0; ii < nin; ++ii) {
      
      // squared distance
      dd     = pow(x - input_xyv(ii, 0), 2); 
      dd    += pow(y - input_xyv(ii, 1), 2)*asp*asp; 
      
      if (max_distance > 0){  // if negative ignore; if above max_distance -> NA
        if (dd > max_distance) dd = 1.0e100;
      }
      
      update_w_v( w, dd, 
                  v, input_xyv(ii, 2), 
                  &imax, ii, nmean);
                  
    }  // end ii
    
    // weighted averages
    vv(io) = get_value (w, v, nmean);
    
  }  
  
  return vv;
}


// =============================================================================
// =============================================================================
//
// interpolation of time series data
// input is as xy data, 2D time format
// output in xy format
//
// =============================================================================
// =============================================================================

// [[Rcpp::export]]

NumericMatrix interpolate_ts_cpp (
    NumericMatrix input_xy,  //  (nin, 2)
    NumericMatrix input_t,   //  (nt,  nin)
    NumericMatrix output_xy, //  (nout, 2)
    int    nmean,    
    double asp) {

  // ------------------------
  // declaration section
  // ------------------------
  
  double x, y, dd  ;
  
  double  *w = new double[nmean] ;    // weights (based on xy)
  int     *iv = new int[nmean] ;      // position of values in input2D
  
  int  imax = 0;
  int  ii, io, iw, it;           // iterators
  
  int  nin = input_xy.nrow();    // size of input geographic points
  int  nt  = input_t.nrow();     // size of input time points
  int  nout = output_xy.nrow();       
  
  NumericMatrix vv(nout, nt);    // output matrix
  
  // ------------------------
  // calculation 
  // ------------------------
  
  for (io = 0; io < nout; ++io) {  // loop over output
    
    x = output_xy(io, 0);
    y = output_xy(io, 1);
    
    for (ii = 0; ii < nin; ++ii) {
      dd     = pow(x - input_xy(ii, 0), 2); 
      dd    += pow(y - input_xy(ii, 1), 2)*asp*asp; 
      
      if (ii < nmean){
        
        w[ii]  = dd;  // save distance
        iv[ii] = ii;  // and position in input

      } else {   // find maximal value
        
        if (ii == nmean) imax = update_imax(w, nmean);
        
        for (iw = 0; iw < nmean; ++iw) {
          
          if (dd < w[imax]){ // replace
            
            w[imax]  = dd;
            iv[imax] = ii;
            imax = update_imax(w, nmean);
          }
          
        } // iw
        
      } // end else
      
      
    }
    // from distance to scaled inverse distances
    inv_distance (w, nmean);
    
    // weighted average for all time points
    for (it = 0; it < nt; it++){
      dd = 0.;
      for (iw = 0; iw < nmean; ++iw) {      
        dd += w[iw] * input_t(it, iv[iw]);
      }
      vv(io, it) = dd;
    } // it
  }  // end io
  
  // ------------------------
  // results
  // ------------------------
  
  return vv;
  
}


// [[Rcpp::export]]

NumericMatrix interpolate_ts_xy_2D_cpp (
    NumericMatrix input_xy,  //  (nin, 2)
    NumericMatrix input_t,   //  (nt,  nin)
    NumericVector output_x,  //  (nox)
    NumericVector output_y,  //  (noy)
    int    nmean,    
    double asp) {
  
  // ------------------------
  // declaration section
  // ------------------------
  
  double x, y, dd  ;
  
  double  *w  = new double[nmean] ;   // weights (based on xy)
  int     *iv = new int[nmean] ;      // position of values in input2D
  
  int  iin, iox, ioy, iw, it, io;         // iterators
  
  int  nin  = input_xy.nrow();    // size of input geographic points
  int  nt   =  input_t.nrow();    // size of input time points
  int  nox  = output_x.size();       
  int  noy  = output_y.size();       
  int  nout = nox*noy;  
  int  imax = 0;
  
  double  *D_x = new double[nin];
  double  D_y;
  
  NumericMatrix vv(nout, nt);    // output matrix
  
  // ------------------------
  // calculation 
  // ------------------------
  io = 0;
  for (iox = 0; iox < nox; ++iox) {  // output x
    
    // squared distance in x for all iin
    x = output_x(iox);
    
    for (iin = 0; iin < nin; ++iin) {  // input
      D_x[iin] = pow(input_xy(iin, 0) - x, 2); 
    }  
    
    for (ioy = 0; ioy < noy; ++ioy) {  // output y
      
      for (iin = 0; iin < nin; ++iin) {  // input
        
        y = input_xy(iin, 1);
        D_y = pow(y - output_y(ioy), 2)*asp*asp; 
        
        dd = D_x[iin] + D_y;  // squared total distance
        
        if (iin < nmean){
          
          w[iin]  = dd;   // save distance
          iv[iin] = iin;  // and position in input
          
        } else {   // find maximal value
          
          if (iin == nmean) imax = update_imax(w, nmean);
          
          for (iw = 0; iw < nmean; ++iw) {
            
            if (dd < w[imax]){ // replace
              
              w[imax]  = dd;
              iv[imax] = iin;
              imax = update_imax(w, nmean);
            }
            
          } // iw
          
        } // end else
      }
      // from distance to scaled inverse distances
      inv_distance (w, nmean);
      
      // weighted average for all time points
      for (it = 0; it < nt; it++){
        dd = 0.;
        for (iw = 0; iw < nmean; ++iw) {      
          dd += w[iw] * input_t(it, iv[iw]);
        }
        vv(io, it) = dd;
      } // it
      io += 1;
    }  // end ioy
  }  // end iox
  
  // ------------------------
  // results
  // ------------------------
  
  return vv;
  
}

