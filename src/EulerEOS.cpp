#include "EulerEOS.h"

EulerEOS::EulerEOS(double gamma)
  :gamma(gamma){};

Eigen::ArrayXXf EulerEOS::prim_to_con(Eigen::ArrayXXf u_p)
{
  Eigen::ArrayXXf u_c(u_p.rows(), u_p.cols());

  u_c.col(0) = u_p.col(0);
  u_c.col(1) = u_p.col(0) * u_p.col(1);
  u_c.col(2) = u_p.col(2) / (gamma - 1.0) + 0.5 * u_p.col(0) * pow(u_p.col(1), 2);

  return u_c;
}

Eigen::ArrayXXf EulerEOS::con_to_prim(Eigen::ArrayXXf u_c)
{
  Eigen::ArrayXXf u_p(u_c.rows(), u_c.cols());

  u_p.col(0) = u_c.col(0);
  u_p.col(1) = u_c.col(1) / u_c.col(0);
  u_p.col(2) = (u_c.col(2) - 0.5 * u_c.col(0) * pow(u_p.col(1), 2.0)) * (gamma - 1.0);
  
  return u_p;
}
Eigen::ArrayXXf EulerEOS::Euler_flux_fn(Eigen::ArrayXXf f, Eigen::ArrayXXf f_prim)
{
  
  Eigen::ArrayXXf flux_fn(f.rows(), f.cols());
  
  for (int i=0; i<f.rows(); i++)
    {
      double rho = f(i, 0);
      double energy = f(i, 2);
      double velocity = f_prim(i, 1);
      double pressure = f_prim(i, 2);
      
      flux_fn(i, 0) = flux_fn_rho(rho, velocity);
      flux_fn(i, 1) = flux_fn_mom(rho, velocity, pressure);
      flux_fn(i, 2) = flux_fn_E(energy, velocity, pressure);
    }

  return flux_fn;

}

double EulerEOS::flux_fn_rho(double rho, double v)
{
  // flux function for density
  return rho * v;
}

double EulerEOS::flux_fn_mom(double rho, double v, double p)
{
  // flux function for momentum
  return rho * pow(v, 2.0) + p;
}

double EulerEOS::flux_fn_E(double E, double v, double p)
{
  // flux function for total energy
  return (E + p) * v;
}
