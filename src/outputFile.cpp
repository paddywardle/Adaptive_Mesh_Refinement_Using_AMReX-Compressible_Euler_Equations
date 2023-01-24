#include "header.H"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>

void Euler1D::outputFile(std::string outputName)
{
  std::ofstream output(outputName);

  for (int i=0; i<u_prim.rows(); i++)
    {
      double x = x0 + (i-1) * dx;
      output<<x<<" ";
      for (int j=0; j<u_prim.cols(); j++)
	{
	  output<<u_prim(i, j)<<" ";
	}
      output<<std::endl;
    }
}
