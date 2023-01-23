#include "header.H"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>

double Euler1D::flux_fn_rho(double rho, double v)
{
  return rho * v;
}

double Euler1D::flux_fn_mom(double rho, double v, double p)
{
  return rho * pow(v, 2.0) + p;
}

double Euler1D::flux_fn_E(double E, double v, double p)
{
  return (E + p) * v;
}

double Euler1D::half_time_step_flux(std::array<double, 3> u_i, std::array<double, 3>u_prim_i, int col)
{
  double f_u;
  
  if (col == 0)
    {
      f_u = flux_fn_rho(u_i[0], u_prim_i[1]);
    }
  else if (col == 1)
    {
      f_u = flux_fn_mom(u_i[0], u_prim_i[1], u_prim_i[2]);
    }
  else if (col == 2)
    {
      f_u = flux_fn_E(u_i[2], u_prim_i[1], u_prim_i[2]);
    }

  return f_u;

}

double Euler1D::lax_friedrich_flux(std::array<double, 3> u_i, std::array<double, 3>u_iPlus1, std::array<double, 3>u_prim_i, std::array<double, 3> u_prim_iPlus1, int col, double dt)
{
  double fhalf;
  
  if (col == 0)
    {
      fhalf = 0.5 * (dx / dt) * (u_i[0] - u_iPlus1[0]) + 0.5 * (flux_fn_rho(u_iPlus1[0], u_prim_iPlus1[1]) + flux_fn_rho(u_i[0], u_prim_i[1]));
    }
  else if (col == 1)
    {
      fhalf = 0.5 * (dx / dt) * (u_i[1] - u_iPlus1[1]) + 0.5 * (flux_fn_mom(u_iPlus1[0], u_prim_iPlus1[1], u_prim_iPlus1[2]) + flux_fn_mom(u_i[0], u_prim_i[1], u_prim_i[2]));
    }
  else if (col == 2)
    {
      fhalf = 0.5 * (dx / dt) * (u_i[2] - u_iPlus1[2]) + 0.5 * (flux_fn_E(u_iPlus1[2], u_prim_iPlus1[1], u_prim_iPlus1[2]) + flux_fn_E(u_i[2], u_prim_i[1], u_prim_i[2]));
    }

  return fhalf;
}

double Euler1D::richtmyer_flux(std::array<double, 3> u_i, std::array<double, 3>u_iPlus1, std::array<double, 3>u_prim_i, std::array<double, 3> u_prim_iPlus1, int col, double dt)
{
  double uhalf_rho = 0.5 * (u_i[0] + u_iPlus1[0]) - 0.5 * (dt / dx) * (flux_fn_rho(u_iPlus1[0], u_prim_iPlus1[1]) - flux_fn_rho(u_i[0], u_prim_i[1]));

  double uhalf_mom = 0.5 * (u_i[1] + u_iPlus1[1]) - 0.5 * (dt / dx) * (flux_fn_mom(u_iPlus1[0], u_prim_iPlus1[1], u_prim_iPlus1[2]) - flux_fn_mom(u_i[0], u_prim_i[1], u_prim_i[2]));

  double uhalf_E = 0.5 * (u_i[2] + u_iPlus1[2]) - 0.5 * (dt / dx) * (flux_fn_E(u_iPlus1[2], u_prim_iPlus1[1], u_prim_iPlus1[2]) - flux_fn_E(u_i[2], u_prim_i[1], u_prim_i[2]));

  double uhalf_v = uhalf_mom / uhalf_rho;

  double uhalf_p = (uhalf_E - 0.5 * uhalf_rho * pow(uhalf_v, 2.0)) * (gamma - 1);

  double fhalf;
  
  if (col == 0)
    {
      fhalf = flux_fn_rho(uhalf_rho, uhalf_v);
    }
  else if (col == 1)
    {
      fhalf = flux_fn_mom(uhalf_rho, uhalf_v, uhalf_p);
    }
  else if (col == 2)
    {
      fhalf = flux_fn_E(uhalf_E, uhalf_v, uhalf_p);
    }

  return fhalf;
}

double Euler1D::FORCE_flux(std::array<double, 3> u_i, std::array<double, 3>u_iPlus1, std::array<double, 3>u_prim_i, std::array<double, 3> u_prim_iPlus1, int col, double dt)
{
  double fhalf = 0.5 * (lax_friedrich_flux(u_i, u_iPlus1, u_prim_i, u_prim_iPlus1, col, dt) + richtmyer_flux(u_i, u_iPlus1, u_prim_i, u_prim_iPlus1, col, dt));

  return fhalf;
}

double Euler1D::calculate_timestep()
{
  
  std::vector<double> wave_speed;
  wave_speed.resize(nCells+2);

  u_prim = con_to_prim(u);
  
  for (int i=0; i<u.size(); i++)
    {
      double cs = sqrt((gamma*u_prim[i][2])/u_prim[i][0]);
      wave_speed[i] = abs(u_prim[i][1]) + cs;
    }

  double a_max = *std::max_element(wave_speed.begin(), wave_speed.end());

  double dt = C * (dx / a_max);

  return dt;
}

double Euler1D::deltai_func(double u_i, double u_iPlus1, double u_iMinus1, double w=0)
{
  return 0.5 * (1 + w) * deltai_half(u_iMinus1, u_i); + 0.5 * (1 - w) * deltai_half(u_i, u_iPlus1);
}

double Euler1D::deltai_half(double u_i, double u_iPlus1)
{
  return u_iPlus1 - u_i;
}

double Euler1D::reconstruction_uL(double u_i, double Xi, double deltai)
{
  return u_i - 0.5 * Xi * deltai;
}

double Euler1D::reconstruction_uR(double u_i, double Xi, double deltai)
{
  return u_i + 0.5 * Xi * deltai;
}

double Euler1D::slope_limiter(double u_i, double u_iPlus1, double u_iMinus1)
{
  //std::cout<<deltai_func(u_iMinus1, u_i) / deltai_func(u_i, u_iPlus1)<<std::endl;
  //std::cout<<u_i<<" "<<u_iPlus1<<" "<<u_iMinus1<<std::endl;
  return deltai_half(u_iMinus1, u_i) / deltai_half(u_i, u_iPlus1);
}

double Euler1D::reconstruction_XiL(double r)
{
  return (2.0 * r) / (1 + r);
}

double Euler1D::reconstruction_XiR(double r)
{
  return 2.0 / (1 + r);
}

double Euler1D::minibee(double r, double XiR)
{
  if (r <= 0)
    {
      return 0.0;
    }
  else if (r > 1)
    {
      return std::min(1.0, XiR);
    }
  else
    {
      return r;
    }
}

void Euler1D::solvers()
{
  initial_conds();
  
  double t = tStart;
  
  std::vector<std::array<double, 3>> flux;
  flux.resize(nCells+1);

  std::vector<std::array<double, 3>> r;
  r.resize(nCells);

  std::vector<std::array<double, 3>> XiL;
  XiL.resize(nCells);

  std::vector<std::array<double, 3>> XiR;
  XiR.resize(nCells);

  std::vector<std::array<double, 3>> Xi;
  Xi.resize(nCells);

  std::vector<std::array<double, 3>> uL;
  uL.resize(nCells);

  std::vector<std::array<double, 3>> uR;
  uR.resize(nCells);

  std::vector<std::array<double, 3>> uL_prim;
  uL_prim.resize(nCells);

  std::vector<std::array<double, 3>> uR_prim;
  uR_prim.resize(nCells);

  std::vector<std::array<double, 3>> uLhalf;
  uLhalf.resize(nCells);

  std::vector<std::array<double, 3>> uRhalf;
  uRhalf.resize(nCells);

  std::vector<std::array<double, 3>> uLhalf_prim;
  uLhalf_prim.resize(nCells);

  std::vector<std::array<double, 3>> uRhalf_prim;
  uRhalf_prim.resize(nCells);
  
  std::vector<std::array<double, 3>> deltai;
  deltai.resize(nCells);
  
  do {
    
    for (int i=1; i<nCells; i++)
      {
	for (int j=0; j<r[i].size(); j++)
	  {
	    r[i-1][j] = slope_limiter(u[i][j], u[i+1][j], u[i-1][j]);
	    XiL[i-1][j] = reconstruction_XiL(r[i-1][j]);
	    XiR[i-1][j] = reconstruction_XiR(r[i-1][j]);
	    Xi[i-1][j] = minibee(r[i-1][j], XiR[i-1][j]);
	    deltai[i-1][j] = deltai_func(u[i][j], u[i+1][j], u[i-1][j]);
	  }
      }
    
    double dt = calculate_timestep();

    u_prim = con_to_prim(u);

    t += dt;

    // add transmissive boundary
    for (int i=0; i<u[0].size(); i++)
      {
	u[0][i] = u[1][i];
	u[nCells+1][i] = u[nCells][i];
      }

    for (int i=0; i<uL.size(); i++)
      {
	for (int j=0; j<uL[i].size(); j++)
	  {
	    uL[i][j] = reconstruction_uL(u[i][j], Xi[i][j], deltai[i][j]);
	    uR[i][j] = reconstruction_uR(u[i][j], Xi[i][j], deltai[i][j]);
	  }
      }

    uL_prim = con_to_prim(uL);
    uR_prim = con_to_prim(uR);

    for (int i=0; i<uL.size(); i++)
      {
	for (int j=0; j<uL[i].size(); j++)
	  {
	    uLhalf[i][j] = uL[i][j] - 0.5 * (dt / dx) * (half_time_step_flux(uR[i], uR_prim[i], j) - half_time_step_flux(uL[i], uL_prim[i], j));
	    uRhalf[i][j] = uR[i][j] - 0.5 * (dt / dx) * (half_time_step_flux(uR[i], uR_prim[i], j) - half_time_step_flux(uL[i], uL_prim[i], j));
	  }
      }

    uLhalf_prim = con_to_prim(uLhalf);
    uRhalf_prim = con_to_prim(uRhalf);
    
    for (int i=0; i<flux.size(); i++)
      {
	for (int j=0; j<flux[i].size(); j++)
	  {
	    flux[i][j] = FORCE_flux(uLhalf[i], uRhalf[i], uLhalf_prim[i], uRhalf_prim[i], j, dt);
	  }
      }
    
    for (int i = 1; i<nCells+1; i++)
      {
	for (int j=0; j<u[i].size(); j++)
	  {
	    uPlus1[i][j] = u[i][j] - (dt/dx) * (flux[i][j] - flux[i-1][j]);
	  }
      }
    u = uPlus1;
    break;
  } while (t < tStop);
}
