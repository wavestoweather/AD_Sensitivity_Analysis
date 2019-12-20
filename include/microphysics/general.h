#ifndef GENERAL_H
#define GENERAL_H

#include <stdio.h>
#include "codi.hpp"

////////////////////////////////////////////////////////////////////////////////
// =============================================================================
// This file contains several general purpose functions.
// =============================================================================
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// This function prints the contents of a vector y with length
// num_comp to stdout.
// 
void printy(const double y[])
{
  std::cout << "\n";
  for(int ii = 0 ; ii < num_comp ; ii++){
    std::cout << "y[" << std::to_string(ii) << "] = " << y[ii] << "\n";
  }
  std::cout << std::endl;
} // End of printy
////////////////////////////////////////////////////////////////////////////////









////////////////////////////////////////////////////////////////////////////////
// This function computes
// y = v1 + a*v2
// for vectors y, v1, v2 and a scalar a.
// 
void v1pav2(double y[],
	    const double v1[],
	    const double v2[],
	    const double a)
{

  for(int ii = 0 ; ii < num_comp ; ii++){
    y[ii] = v1[ii] + a*v2[ii];
  }

} // End of v1pav2
////////////////////////////////////////////////////////////////////////////////

#endif
