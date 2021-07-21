#include <cosDispersionC.h>
#include <dispersionFitC.h>
#include <stdio.h>
#include <math.h>

int main() {
  float x;
  // call a function in another file
  x = linear_fit(10, 1, 2);
  printf("%.6f", x);
  return(0);
}