#include "header.H"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>

std::vector<std::array<double, 3>> Euler1D::prim_to_con(std::vector<std::array<double, 3>> u_p)
{
  std::vector<std::array<double, 3>> u_c;
  u_c.resize(nCells+2);
  
  for (int i=0; i<u_p.size(); i++)
    {
      u_c[i][0] = u_p[i][0];
      u_c[i][1] = u_p[i][0] * u_p[i][1];
      u_c[i][2] = u_p[i][2]/(gamma-1.0) + 0.5 * u_p[i][0] * pow(u_p[i][1], 2.0);
    }

  return u_c;
}

std::vector<std::array<double, 3>> Euler1D::con_to_prim(std::vector<std::array<double, 3>> u_c)
{
  std::vector<std::array<double, 3>> u_p;
  u_p.resize(nCells+2);
  
  for (int i=0; i<u_c.size(); i++)
    {
      u_p[i][0] = u_c[i][0];
      u_p[i][1] = u_c[i][1] / u_c[i][0];
      u_p[i][2] = (u_c[i][2] - 0.5 * u_p[i][0]*pow(u_p[i][1],2.0)) * (gamma - 1.0);
    }

  return u_p;
}
