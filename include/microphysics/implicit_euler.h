#ifndef IMPLICIT_EULER_H
#define IMPLICIT_EULER_H

#include "codi.hpp"

#include "constants.h"
#include "user_functions.h"



////////////////////////////////////////////////////////////////////////////////
// =============================================================================
// This file provides the function for a single step with the
// Implicit Euler method
// =============================================================================
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
// This method implements the computation of the quadratic norm of a
// difference
// y1-y1
// for two given vectors y1 and y2.
//
inline codi::RealReverse compute_diff_length_square( codi::RealReverse y1[],
						     codi::RealReverse y2[],
						     const int num_comp )
{
  codi::RealReverse result = 0.0;
  
  for(int ii = 0 ; ii < num_comp ; ii++){
    result = result + (y1[ii]-y2[ii])*(y1[ii]-y2[ii]);
  }

  return result;
} // End of method "compute_diff_length_square"
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
// This method computes a single step using the implicit method of
// Euler for the ODE
// y' = RHS(y)
// with old system state yold and new
// system state ynew.  Using the implicit method of Euler involves the
// solution of an implicit equation in each step. In this
// implementation, we simply solve this implicit equation using a
// fixed-point iteration. This has the drawback, that the method
// should _NOT_ be used for stiff equations.
//
void implicit_euler_step( codi::RealReverse ynew[],
			  codi::RealReverse yold[],
			  const int num_comp,
			  const reference_quantities_t& ref,
			  model_constants_t& cc)
{

  // Define the tolerance
  const double TOL_square = (1.0e-8)*(1.0e-8);
  codi::RealReverse err_square = 1.0e10;

  // Define the temporary variables
  codi::RealReverse yit[num_comp];
  codi::RealReverse k[num_comp];
  int count = 0;


  // Evaluate the RHS to get f(yold)
  RHS(k, yold, num_comp, ref, cc);

  // Compute the first guess
  for(int ii = 0 ; ii < num_comp ; ii++){
    yit[ii] = yold[ii] + cc.dt_half*k[ii];
  }


  while( err_square > TOL_square )
    {
      // Safety measure
      if(count > 10000){
	std::cout << "\nWARNING: Reached the maximal step count for the fixed-point iteration!\n";
	break;
      }
      count++;

      // Compute the new iterate
      RHS(k, yit, num_comp, ref, cc);
      for(int ii = 0 ; ii < num_comp ; ii++){
	ynew[ii] = yold[ii] + cc.dt_half*k[ii];
      }

      // Compute the error
      err_square = compute_diff_length_square(ynew, yit, num_comp);

      // 
      // Swap the contents of vectors yit and ynew
      // NOTE: This operation could be done using pointers. BUT: is it also possible with AD?
      //
      for(int ii = 0 ; ii < num_comp ; ii++){
	yit[ii] = ynew[ii];
      }
      
      
    }
  

} // End of method implicit_euler_step
////////////////////////////////////////////////////////////////////////////////









#endif
