#include <Rcpp.h>
using namespace Rcpp;

// 
// C++ functions to average dense points into a coarser grid
//

// [[Rcpp::plugins("cpp11")]]

// =============================================================================
// =============================================================================
// Maps xyt data from a very high resolution, unstructured grid to a 
// lower resolution regular grid by averaging in each grid cell.
// =============================================================================
// =============================================================================

// [[Rcpp::export]]

NumericMatrix average_grid_t(      // averaged grid
    
    NumericMatrix   in_xy2D, // input matrix: x, y, multiple values, one for each time  
    NumericVector   in_t,    // input times - corresponds to in_xy2D (columns 3- last..)
    NumericMatrix   xyout,   // xy output values
    NumericVector   range_x, // x grid - 3 values: min x, max x and dx
    NumericVector   range_y, // y grid - 3 values: min y, max y and dy
    NumericVector   out_t)   // output times - vector
  
{  
  
  // ========================
  // declaration section
  // ========================
  
  int  nr = in_xy2D.nrow();      // number of input data rows
  int  nt = in_xy2D.ncol() - 2;  // number of input time points (first two columns: x and y)
  
  if (in_t.size() != nt)
    stop("Error: in_t should be as long as the number of columns of in_xy2D - 2");
  
  // int nt = in_t.size();       // should be the same
  int  nout = xyout.nrow();      // number of output data rows
  
  int  i, ir, it, ii, io;        // iterators
  int  ix, iy;                   // index to x and y of input data
  
  double tt, dt, tfac, x, y ;
  bool hasData;
  
  // output grid in x
  double xmin = range_x(0);
  double xmax = range_x(1);
  double dx   = range_x(2);
  
  int nx   = int((xmax - xmin)/dx) + 1;
  
  // output grid in y
  double ymin = range_y(0);
  double ymax = range_y(1);
  double dy   = range_y(2);
  
  int ny   = int((ymax - ymin)/dy) + 1;
  
  // 2D matrix with output values
  IntegerMatrix i2D(nx, ny);
  
  for (ix = 0; ix < nx; ix++){
    for (iy = 0; iy < ny; iy++){
      i2D(ix, iy) = -99;
    }
  }
    
  // position of xyout x and y in 2D output grid. Add the number to the grid cell 
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
  
  // ========================
  // loop over all input points
  // ========================
  
  // for each input: find grid cell, add values; count the number of values 
  
  for (ir = 0; ir < nr; ++ir) {
    
    // find ix, iy 
    x = in_xy2D(ir, 0);
    y = in_xy2D(ir, 1);
    
    ix = int((x-xmin)/dx);  // index to 2D matrix in x
    iy = int((y-ymin)/dy);  //                   and y
    
    if ((ix < nx) & (ix >= 0) & 
        (iy < ny) & (iy >= 0)){

      ii = i2D(ix, iy);
      
      if ( ii >= 0) { // it is one of the requested output xy
        
        for (it = 0; it < nt; ++it){
          average  (ii, it) += in_xy2D(ir, 2+it);   // add value
          numpoints(ii, it) += 1;                   // add 1
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
    
    tt = in_t(it);
    dt = in_t(it+1) - tt;

    while (out_t(io) <= in_t(it+1)) {

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


