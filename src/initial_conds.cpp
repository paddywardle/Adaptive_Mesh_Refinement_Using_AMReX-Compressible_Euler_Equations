#include "header.H"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>

void Euler1D::initial_conds()
{
  for (int i=0; i<u.rows(); i++)
    {
      double x = x0 + i * dx;
      
      if (x <= 0.5)
	{
	  // density
	  u_prim(i, 0) = 1.0;
	  // velocity
	  u_prim(i, 1) = 0.0;
	  // pressure
	  u_prim(i, 2) = 1.0;
	}
      else
	{
	  // density
	  u_prim(i, 0) = 0.125;
	  // velocity
	  u_prim(i, 1) = 0.0;
	  // pressure
	  u_prim(i, 2) = 0.1;
	}
    }
}
