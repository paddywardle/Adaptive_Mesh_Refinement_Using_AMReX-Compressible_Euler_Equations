#include "header.H"

void Euler1D::outputFile(std::string outputName)
{
  std::ofstream output(outputName);

  for (int i=1; i<u_prim.rows()-2; i++)
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
