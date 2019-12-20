#ifndef EULER_H
#define EULER_H

#include "codi.hpp"
#include "constants.h"
#include "user_functions.h"

////////////////////////////////////////////////////////////////////////////////
// =============================================================================
// This file provides the function for a single step with the
// explicit Euler method
// =============================================================================
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
// This method computes a single Euler-step with timestep dt for the ODE
// y' = RHS(y)
// taking yold as the initial value and storing the result in ynew.
void euler_step(codi::RealReverse ynew[],
		codi::RealReverse yold[],
		const int num_comp,
		const reference_quantities_t &ref, 
		model_constants_t &cc)
{
  codi::RealReverse k1[num_comp];

  // Evaluate the RHS
  RHS(k1, yold, num_comp, ref, cc);

  // Do the step
  for(int ii = 0 ; ii < num_comp ; ii++){
    ynew[ii] = yold[ii] + cc.dt*k1[ii];
  }

  
} // End of method euler_step
////////////////////////////////////////////////////////////////////////////////

#endif
