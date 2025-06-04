#include <Rcpp.h>
using namespace Rcpp;

// 
// C++ functions to average dense points into a coarser grid
//

// [[Rcpp::plugins("cpp11")]]

// =============================================================================
// =============================================================================
// Maps xy data from a very high resolution, unstructured grid to a 
// lower resolution regular grid by averaging in each grid cell.
// only for regular output grids.
// =============================================================================
// =============================================================================

// [[Rcpp::export]]

NumericMatrix average_xy_2D_r_cpp( 
      NumericMatrix  input_xyv,
      NumericVector range_x, 
      NumericVector range_y){
  // ------------------------
  // declaration section
  // ------------------------
  
  double xi, yi, v  ;

  int  iox, ioy, iin;                   // iterators
  
  // output grid in x
  double xmin = range_x(0);  
  double xmax = range_x(1);  
  double dx   = range_x(2);
  
  // output grid in y
  double ymin = range_y(0);
  double ymax = range_y(1);  
  double dy   = range_y(2);
  
  int  nin   = input_xyv.nrow();       

  int  nox =  int ((xmax - xmin)/dx + 0.01) +1;       
  int  noy =  int ((ymax - ymin)/dy + 0.01) +1;       
  
  IntegerMatrix   indx(nox, noy);      // counts
  NumericMatrix   vv  (nox, noy);      // output values
  

  // ------------------------
  // calculation 
  // ------------------------
  
  // initialisation
  
  for (iox = 0; iox < nox; ++iox) {
    for (ioy = 0; ioy < noy; ++ioy) {
      indx(iox, ioy) = 0;
      vv(iox, ioy) = 0.;
    }
  }
  
  // sum the inputs in the output grid
  
  for (iin = 0; iin < nin; ++iin) {
    
    xi = input_xyv(iin, 0);
    yi = input_xyv(iin, 1);
    v  = input_xyv(iin, 2);
      
    if ((v != NA_REAL) &
        (xi >= xmin) & (xi <= xmax) & 
        (yi >= ymin) & (yi <= ymax)){
    
    // position in output grid, x-direction
       int ix = (xi - xmin)/dx;
       int iy = (yi - ymin)/dy;
       
       indx(ix, iy) ++ ;
       vv  (ix, iy) += input_xyv(iin, 2) ;
       
    }
  }
  
  // calculate the mean
  
  for (iox = 0; iox < nox; ++iox) {
    for (ioy = 0; ioy < noy; ++ioy) {
      if (indx(iox, ioy) > 0){
        vv(iox, ioy) = vv(iox, ioy)/indx(iox, ioy);    
      }  else {
        vv(iox, ioy) = NA_REAL;    
         }
      }
    }  
  return vv;  
} 

// =============================================================================
// =============================================================================
// Maps xy data from a very high resolution, regular grid to a 
// lower resolution regular grid by averaging in each grid cell.
// only for regular output grids.
// =============================================================================
// =============================================================================

// [[Rcpp::export]]

NumericMatrix average_2D_2D_r_cpp(
      NumericVector   input_x,
      NumericVector   input_y,
      NumericMatrix  input_2D,
      NumericVector range_x, 
      NumericVector range_y){

    double xi, yi, v  ;
    
    int  iox, ioy, iix, iiy;                   // iterators
    
    // output grid in x
    double xmin = range_x(0);  
    double xmax = range_x(1);  
    double dx   = range_x(2);
    
    // output grid in y
    double ymin = range_y(0);
    double ymax = range_y(1);  
    double dy   = range_y(2);
    
    int  nix   = input_2D.nrow();       
    int  niy   = input_2D.ncol();       
    
    int  nox =  int ((xmax - xmin)/dx + 0.01) +1;       
    int  noy =  int ((ymax - ymin)/dy + 0.01) +1;       
    
    IntegerMatrix   indx(nox, noy);      // counts
    NumericMatrix   vv  (nox, noy);      // output values
    
    
    // ------------------------
    // calculation 
    // ------------------------
    
    // initialisation
    
    for (iox = 0; iox < nox; ++iox) {
      for (ioy = 0; ioy < noy; ++ioy) {
        indx(iox, ioy) = 0;
        vv(iox, ioy) = 0.;
      }
    }
    
    // sum the inputs in the output grid
    
    for (iix = 0; iix < nix; ++iix) {
      
      xi = input_x(iix);
      
      if ((xi >= xmin) & (xi <= xmax)){ 
        
        for (iiy = 0; iiy < niy; ++iiy) {
          yi = input_y(iiy);
          v  = input_2D(iix, iiy);  
          if ((v != NA_REAL) &
              (yi >= ymin) & (yi <= ymax)){
        
           // position in output grid, x-direction
            int ix = (xi - xmin)/dx;
            int iy = (yi - ymin)/dy;
        
            indx(ix, iy) ++  ;
            vv  (ix, iy) += v;
          }
        }
      }
   }    
    
    // calculate the mean
    for (iox = 0; iox < nox; ++iox) {
      for (ioy = 0; ioy < noy; ++ioy) {
        if (indx(iox, ioy) > 0){
          vv(iox, ioy) = vv(iox, ioy)/indx(iox, ioy);    
        }  else {
          vv(iox, ioy) = NA_REAL;    
        }
      }
    }  
    return vv;  
} 

// =============================================================================
// =============================================================================
// **** TIME-SERIES ****
// Maps xyt data from a very high resolution, unstructured grid to a 
// lower resolution regular grid by averaging in each grid cell.
// =============================================================================
// =============================================================================

// [[Rcpp::export]]

NumericMatrix average_grid_t(      // averaged grid
    
    NumericMatrix   input_xy,   // input matrix: x, y  
    NumericMatrix   input_2Dt,  // input matrix: time-series values  
    NumericVector   input_t,    // input times - corresponds to input_xy2D (columns 3- last..)
    NumericMatrix   xyout,   // xy output values
    NumericVector   outx, // output x grid - 3 values: min x, max x and dx
    NumericVector   outy, // output y grid - 3 values: min y, max y and dy
    NumericVector   out_t)   // output times - vector
  
{  
  
  // ========================
  // declaration section
  // ========================
  
  int  nr = input_xy.nrow();     // number of input data rows (x, y) values
  int  nt = input_2Dt.ncol();    // number of input time points 
  
  if (input_t.size() != nt)
    stop("Error: input_t should be as long as the number of columns of input_2Dt");
  
  if (input_2Dt.nrow() != nr)
    stop("Error: input_2Dt should have as many rows as input_xy");
  
  // int nt = input_t.size();       // should be the same
  int  nout = xyout.nrow();      // number of output data rows
  
  int  i, ir, it, ii, io;        // iterators
  int  ix, iy;                   // index to x and y of input data
  
  double tt, dt, tfac, x, y ;
  bool hasData;
  
  // output grid in x
  double xmin = outx(0);
  double xmax = outx(1);
  double dx   = outx(2);
  
  int nx   = int((xmax - xmin)/dx) + 1;
  
  // output grid in y
  double ymin = outy(0);
  double ymax = outy(1);
  double dy   = outy(2);
  
  int ny   = int((ymax - ymin)/dy) + 1;
  
  // 2D matrix with output values
  IntegerMatrix i2D(nx, ny);
  
  
  // =========================================================
  // Find if there is request for output in a grid cell, based on xyout
  // =========================================================
  
  // initialise
  for (ix = 0; ix < nx; ix++){
    for (iy = 0; iy < ny; iy++){
      i2D(ix, iy) = -99;
    }
  }
  
  // loop over all nout cells
  for (i = 0; i < nout; i++){
    x = xyout(i, 0) + dx/100.;   // add dx/100 to avoid roundoff
    y = xyout(i, 1) + dy/100.;
    
    ix = int((x-xmin)/dx);
    iy = int((y-ymin)/dy);
    
    if ((ix < nx) & (ix >= 0) & 
        (iy < ny) & (iy >= 0)){
      i2D(ix, iy) = i;
    }  
  }
  
  NumericMatrix average   (nout, nt);
  IntegerMatrix numpoints (nout, nt); 
  
  // =========================================================
  // loop over all input points and create input time-series
  // =========================================================
  
  // for each input: find grid cell, add values; count the number of values 
  
  for (ir = 0; ir < nr; ++ir) {
    
    // find ix, iy 
    x = input_xy(ir, 0);
    y = input_xy(ir, 1);
    
    ix = int((x-xmin)/dx);  // index to 2D matrix in x
    iy = int((y-ymin)/dy);  //                   and y
    
    if ((ix < nx) & (ix >= 0) & 
        (iy < ny) & (iy >= 0)){
      
      ii = i2D(ix, iy);
      
      if ( ii >= 0) { // it is one of the requested output xy
        
        for (it = 0; it < nt; ++it){
          if (input_2Dt(ir, it) != NA_REAL){
            average  (ii, it) += input_2Dt(ir, it);   // add value
            numpoints(ii, it) += 1;                   // add 1
          }
        }            
      }   
    }  
  }
  
  // ========================
  // mean = total / numpoints
  // discard grid cells with 0 values
  // ========================
  
  for (ir = 0; ir < nout; ++ir) {
    
    hasData = true;
    
    for (it = 0; it < nt; ++it){
      if (numpoints(ir, it) == 0) { // one missing value is enough to discard
        hasData = false;
        break;
      }
      average(ir, it) = average(ir, it)/ numpoints(ir, it);   
    }
    
    if (!hasData) {
      for (it = 0; it < nt; ++it) {
        average(ir, it) = -99.;
      }
    }  
  } 
  
  // =========================================================
  // loop over all output points and create output time-series
  // =========================================================
  
  int ntout  = out_t.size();  // number of output time points
  
  NumericMatrix out  (nout, 2 + ntout); 
  
  // =================================================
  // save the x and y positions in first two columns
  // =================================================
  
  for (ir = 0; ir < nout; ++ir) {
    out(ir, 0) = xyout(ir, 0);
    out(ir, 1) = xyout(ir, 1);
  }  
  
  // ========================
  // approx to output time
  // ========================
  
  io = 0;
  it = 0;
  
  for (it = 0; it < nt - 1; ++it){  
    
    if (io >= ntout) break;
    
    tt = input_t(it);
    dt = input_t(it+1) - tt;
    
    while (out_t(io) <= input_t(it+1)) {
      
      tfac = (out_t(io) - tt)/dt;
      
      for (ii = 0; ii < nout; ++ii){
        out(ii, 2 + io) = average(ii, it) + tfac*( average(ii, it+1)-  average(ii, it));
      }
      io++;
      if (io >= ntout) break;
      
    }   
  }
  
  return out;
}
