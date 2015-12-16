#include "math.h"
// CDF poisson
// compute the probability of P(y < x) of Poisson distribution with mu.
double cdf_poisson(int x, float mu) {
   int i;
   double px, cx;
   px = 1.0;
   cx = px;
   for(i=1; i<=x; i++) {
     px = px * mu / i;
     cx += px;
   }
   return cx* exp(-mu);

}

