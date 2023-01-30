#include "header.H"

Eigen::ArrayXXf Euler1D::prim_to_con(Eigen::ArrayXXf u_p)
{
  Eigen::ArrayXXf u_c(u_p.rows(), u_p.cols());

  u_c.col(0) = u_p.col(0);
  u_c.col(1) = u_p.col(0) * u_p.col(1);
  u_c.col(2) = u_p.col(2) / (gamma - 1) + 0.5 * u_p.col(0) * pow(u_p.col(1), 2);

  return u_c;
}

Eigen::ArrayXXf Euler1D::con_to_prim(Eigen::ArrayXXf u_c)
{
  Eigen::ArrayXXf u_p(u_c.rows(), u_c.cols());

  u_p.col(0) = u_c.col(0);
  u_p.col(1) = u_c.col(1) / u_c.col(0);
  u_p.col(2) = (u_c.col(2) - 0.5 * u_c.col(0) * pow(u_p.col(1), 2.0)) * (gamma - 1.0);
  
  return u_p;
}
