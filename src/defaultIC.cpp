#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector defaultIC(NumericMatrix ppm, NumericVector bkg) {
   int nrow = ppm.nrow(), ncol = ppm.ncol();
   NumericMatrix ppmscore(nrow, ncol)
   ppmscore = ppm *
   NumericVector out(ncol);
   for(int i = 0; i < ncol; i++) {
     double total = 0;
     for(int j = 0; j < nrow; j++) {
       total += ppm(i, j)
     }
   }
}
