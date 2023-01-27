#include "header.H"

Eigen::MatrixXd Euler1D::prim_to_con(Eigen::MatrixXd u_p)
{
  Eigen::MatrixXd u_c(u_p.rows(), u_p.cols());
  
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
  Eigen::MatrixXd u_p(u_c.rows(), u_c.cols());
  /*
  for (int i=0; i<u_c.rows(); i++)
    {
      u_p(i, 0) = u_c(i, 0);
      u_p(i, 1) = u_c(i, 1) / u_c(i, 0);
      u_p(i, 2) = (u_c(i, 2) - 0.5 * u_p(i, 0) * pow(u_p(i, 1),2.0)) * (gamma - 1.0);
    }
  */
  
  u_p.col(0) = u_c.col(0);
  u_p.col(1) = u_c.col(1) / u_c.col(0);
  u_p.col(2) = (u_c.col(2) - 0.5 * u_c.col(0) * pow(u_p.col(1), 2.0)) * (gamma - 1.0);
  
  return u_p;
}
