#include "header.H"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>

Eigen::MatrixXd Euler1D::prim_to_con(Eigen::MatrixXd u_p)
{
  Eigen::MatrixXd u_c(nCells+2, 3);
  
  for (int i=0; i<u_p.rows(); i++)
    {
      u_c(i, 0) = u_p(i, 0);
      u_c(i, 1) = u_p(i, 0) * u_p(i, 1);
      u_c(i, 2) = u_p(i, 2)/(gamma-1.0) + 0.5 * u_p(i, 0) * pow(u_p(i, 1), 2.0);
    }

  return u_c;
}

Eigen::MatrixXd Euler1D::con_to_prim(Eigen::MatrixXd u_c)
{
  Eigen::MatrixXd u_p(nCells+2, 3);
  
  for (int i=0; i<u_c.rows(); i++)
    {
      u_p(i, 0) = u_c(i, 0);
      u_p(i, 1) = u_c(i, 1) / u_c(i, 0);
      u_p(i, 2) = (u_c(i, 2) - 0.5 * u_p(i, 0)*pow(u_p(i, 2),2.0)) * (gamma - 1.0);
    }

  return u_p;
}
