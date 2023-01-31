#include "NumericalMethod.h"

NumericalMethod::NumericalMethod(double gamma, int nCells)
  :gamma(gamma), nCells(nCells), eos(gamma){};

double NumericalMethod::deltai_func(double u_i, double u_iPlus1, double u_iMinus1, double w=0.0)
{
  // calculates cell delta value for boundary extrapolated reconstruction
  return 0.5 * (1.0 + w) * (u_i - u_iMinus1) + 0.5 * (1.0 - w) * (u_iPlus1 - u_i);
}

double NumericalMethod::slope_limiter(double u_i, double u_iPlus1, double u_iMinus1)
{
  // calculates slope ratio value
  
  if ((u_iPlus1 - u_i) == 0.0){
    return 0.0;
  }
  return (u_i - u_iMinus1) / (u_iPlus1 - u_i);
}

double NumericalMethod::limiterXi(double r, Limiters limiter)
{

  double Xi;
  
  if (limiter == Limiters::Superbee)
    {
      Xi = superbee(r);
    }
  else if (limiter == Limiters::Van_Leer)
    {
      Xi = van_leer(r);
    }
  else if (limiter == Limiters::Van_Albada)
    {
      Xi = van_albada(r);
    }
  else
    {
      Xi = minbee(r);
    }

  return Xi;
}

double NumericalMethod::reconstruction_uL(double u_i, double u_iPlus1, double u_iMinus1, Limiters limiter)
{
  double r = slope_limiter(u_i, u_iPlus1, u_iMinus1);
  double Xi = limiterXi(r, limiter);
  double deltai = deltai_func(u_i, u_iPlus1, u_iMinus1);

  return u_i - 0.5 * Xi * deltai;
}

double NumericalMethod::reconstruction_uR(double u_i, double u_iPlus1, double u_iMinus1, Limiters limiter)
{
  double r = slope_limiter(u_i, u_iPlus1, u_iMinus1);
  double Xi = limiterXi(r, limiter);
  double deltai = deltai_func(u_i, u_iPlus1, u_iMinus1);

  return u_i + 0.5 * Xi * deltai;
}

double NumericalMethod::minbee(double r)
{
  if (r <= 0.0)
    {
      return 0.0;
    }
  else if (r > 1.0)
    {
      return std::min(1.0, (2.0 / (1.0 + r)));
    }
  else
    {
      return r;
    }
}

double NumericalMethod::superbee(double r)
{
  if (r <= 0.0)
    {
      return 0.0;
    }
  else if (r > 1.0 && r <= 0.5)
    {
      return 2.0 * r;
    }
  else if (r < 0.5 && r <= 1.0)
    {
      return 1.0;
    }
  else
    {
      return std::min({r, (2.0 / (1.0 + r)), 2.0});
    }
}

double NumericalMethod::van_leer(double r)
{
  if (r <= 0.0)
    {
      return 0.0;
    }
  else
    {
      return std::min({((2.0 * r) / (1.0 + r)), (2.0 / (1.0 + r))});
    }
}

double NumericalMethod::van_albada(double r)
{
  if (r <= 0.0)
    {
      return 0.0;
    }
  else
    {
      return std::min((r*(1.0 + r) / (1 + pow(r, 2.0))), (2.0 / (1.0 + r)));
    }
}

double NumericalMethod::lax_friedrich_flux(double u_i, double u_iPlus1, double flux_i, double flux_Plus1,  double dt, double dx)
{
  
  double fhalf;

  fhalf = 0.5 * (dx / dt) * (u_i - u_iPlus1) + 0.5 * (flux_Plus1 + flux_i);
  
  return fhalf;
}

double NumericalMethod::richtmyer_flux(double u_i, double u_iPlus1, double flux_i, double flux_Plus1,  double dt, double dx)
{

  double uhalf;
  double fhalf;

  uhalf = 0.5 * (u_i + u_iPlus1) - 0.5 * (dt / dx) * (flux_Plus1 - flux_i);
  
  return uhalf;
}

Eigen::ArrayXXf NumericalMethod::FORCE_flux(Eigen::ArrayXXf uLhalf,  Eigen::ArrayXXf uRhalf, Eigen::ArrayXXf uLhalf_prim, Eigen::ArrayXXf uRhalf_prim, double dt, double dx)
{
  
  Eigen::ArrayXXf fhalf(nCells+1, uLhalf.cols());
  Eigen::ArrayXXf uRflux = eos.Euler_flux_fn(uRhalf, uRhalf_prim);
  Eigen::ArrayXXf uLflux = eos.Euler_flux_fn(uLhalf, uLhalf_prim);
  Eigen::ArrayXXf uhalf(nCells+1, uLhalf.cols());
  Eigen::ArrayXXf RI_flux(nCells+1, uLhalf.cols());
  Eigen::ArrayXXf LF_flux(nCells+1, uLhalf.cols());

  for (int i=0; i<nCells+1; i++)
    {
      for (int j=0; j<uLhalf.cols(); j++)
	{
	  uhalf(i, j) = richtmyer_flux(uRhalf(i, j), uLhalf(i+1, j), uRflux(i, j), uLflux(i+1, j), dt, dx);
	  LF_flux(i, j) = lax_friedrich_flux(uRhalf(i, j), uLhalf(i+1, j), uRflux(i, j), uLflux(i+1, j), dt, dx);
	}
    }

  RI_flux = eos.Euler_flux_fn(uhalf, eos.con_to_prim(uhalf));

  for (int i=0; i<nCells+1; i++)
    {
      for (int j=0; j<uLhalf.cols(); j++)
	{
	  fhalf(i, j) = 0.5 * (LF_flux(i, j) + RI_flux(i, j));
	}
    }

  return fhalf;
}

Eigen::ArrayXXf NumericalMethod::uL_half_update(Eigen::ArrayXXf uL, Eigen::ArrayXXf uR, double dt, double dx)
{
  
  // data
  Eigen::ArrayXXf uL_flux = eos.Euler_flux_fn(uL, eos.con_to_prim(uL));
  Eigen::ArrayXXf uR_flux = eos.Euler_flux_fn(uR, eos.con_to_prim(uR));

  return uL - 0.5 * (dt / dx) * (uR_flux - uL_flux);
  
}

Eigen::ArrayXXf NumericalMethod::uR_half_update(Eigen::ArrayXXf uL, Eigen::ArrayXXf uR, double dt, double dx)
{
  
  // data
  Eigen::ArrayXXf uL_flux = eos.Euler_flux_fn(uL, eos.con_to_prim(uL));
  Eigen::ArrayXXf uR_flux = eos.Euler_flux_fn(uR, eos.con_to_prim(uR));

  return uR - 0.5 * (dt / dx) * (uR_flux - uL_flux);
  
}
