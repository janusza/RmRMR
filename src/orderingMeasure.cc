#include "R.h"
#include "Rmath.h"

#include <vector>

extern "C" {
void  computeSwitchesC(double *target, int *N, int *value) {

  for(int i = 0; i < N[0] - 1; ++i)	{
    if(target[i] != target[i+1])	{
      ++value[0];
    }
  }
  return;
}
}
