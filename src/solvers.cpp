#include "header.H"

Eigen::ArrayXXf Euler1D::Euler_flux_fn(Eigen::ArrayXXf f, Eigen::ArrayXXf f_prim)
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

double Euler1D::flux_fn_rho(double rho, double v)
{
  // flux function for density
  return rho * v;
}

double Euler1D::flux_fn_mom(double rho, double v, double p)
{
  // flux function for momentum
  return rho * pow(v, 2.0) + p;
}

double Euler1D::flux_fn_E(double E, double v, double p)
{
  // flux function for total energy
  return (E + p) * v;
}

double Euler1D::calculate_timestep()
{
  
  std::vector<double> wave_speed;
  wave_speed.resize(nCells+4);

  u_prim = con_to_prim(u);
  
  for (int i=0; i<u.rows(); i++)
    {
      double cs = sqrt((gamma*u_prim(i, 2))/u_prim(i, 0));
      wave_speed[i] = abs(u_prim(i, 1)) + cs;
    }

  double a_max = *std::max_element(wave_speed.begin(), wave_speed.end());

  double dt = C * (dx / a_max);

  return dt;
}

double Euler1D::deltai_func(double u_i, double u_iPlus1, double u_iMinus1, double w=0.0)
{
  // calculates cell delta value for boundary extrapolated reconstruction
  return 0.5 * (1.0 + w) * (u_i - u_iMinus1) + 0.5 * (1.0 - w) * (u_iPlus1 - u_i);
}

double Euler1D::reconstruction_uL(double u_i, double Xi, double deltai)
{
  // calculates left boundary reconstruction value
  return u_i - 0.5 * Xi * deltai;
}

double Euler1D::reconstruction_uR(double u_i, double Xi, double deltai)
{
  // calculates right boundary reconstruction value
  return u_i + 0.5 * Xi * deltai;
}

double Euler1D::slope_limiter(double u_i, double u_iPlus1, double u_iMinus1)
{
  // calculates slope ratio value
  
  if ((u_iPlus1 - u_i) == 0.0){
    return 0.0;
  }
  return (u_i - u_iMinus1) / (u_iPlus1 - u_i);
}

double Euler1D::minbee(double r)
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

double Euler1D::lax_friedrich_flux(double u_i, double u_iPlus1, double flux_i, double flux_Plus1,  double dt)
{
  
  double fhalf;

  fhalf = 0.5 * (dx / dt) * (u_i - u_iPlus1) + 0.5 * (flux_Plus1 + flux_i);
  
  return fhalf;
}

double Euler1D::richtmyer_flux(double u_i, double u_iPlus1, double flux_i, double flux_Plus1,  double dt)
{

  double uhalf;
  double fhalf;

  uhalf = 0.5 * (u_i + u_iPlus1) - 0.5 * (dt / dx) * (flux_Plus1 - flux_i);
  
  return uhalf;
}

Eigen::ArrayXXf Euler1D::FORCE_flux(Eigen::ArrayXXf uLhalf,  Eigen::ArrayXXf uRhalf, Eigen::ArrayXXf uLhalf_prim, Eigen::ArrayXXf uRhalf_prim, double dt)
{
  Eigen::ArrayXXf fhalf(nCells+1, uLhalf.cols());
  Eigen::ArrayXXf uRflux = Euler_flux_fn(uRhalf, uRhalf_prim);
  Eigen::ArrayXXf uLflux = Euler_flux_fn(uLhalf, uLhalf_prim);
  Eigen::ArrayXXf uhalf(nCells+1, uLhalf.cols());
  Eigen::ArrayXXf RI_flux(nCells+1, uLhalf.cols());
  Eigen::ArrayXXf LF_flux(nCells+1, uLhalf.cols());

  for (int i=0; i<nCells+1; i++)
    {
      for (int j=0; j<uLhalf.cols(); j++)
	{
	  uhalf(i, j) = richtmyer_flux(uRhalf(i, j), uLhalf(i+1, j), uRflux(i, j), uLflux(i+1, j), dt);
	  LF_flux(i, j) = lax_friedrich_flux(uRhalf(i, j), uLhalf(i+1, j), uRflux(i, j), uLflux(i+1, j), dt);
	}
    }

  RI_flux = Euler_flux_fn(uhalf, con_to_prim(uhalf));

  for (int i=0; i<nCells+1; i++)
    {
      for (int j=0; j<uLhalf.cols(); j++)
	{
	  fhalf(i, j) = 0.5 * (LF_flux(i, j) + RI_flux(i, j));
	}
    }

  return fhalf;
}

void Euler1D::solvers()
{
  // sets initial conditions into instance variables
  initial_conds();

  // set current time to simulation start time
  double t = tStart;

  Eigen::ArrayXXf flux(nCells+1, 3);
  Eigen::ArrayXXf uL(nCells+2, 3);
  Eigen::ArrayXXf uR(nCells+2, 3);
  Eigen::ArrayXXf uL_prim(nCells+2, 3);
  Eigen::ArrayXXf uR_prim(nCells+2, 3);
  Eigen::ArrayXXf uLhalf(nCells+2, 3);
  Eigen::ArrayXXf uRhalf(nCells+2, 3);
  Eigen::ArrayXXf uLhalf_prim(nCells+2, 3);
  Eigen::ArrayXXf uRhalf_prim(nCells+2, 3);

  Eigen::ArrayXXf uL_flux(nCells+2, 3);
  Eigen::ArrayXXf uR_flux(nCells+2, 3);
  
  do {
    
    // calculate timestep based on initial data
    double dt = calculate_timestep();

    // add timestep to current time
    t += dt;
    
    // add transmissive boundary
    for (int i=0; i<u.cols(); i++)
      {
	u(0,i) = u(2, i);
	u(1,i) = u(2, i);
	u(nCells+2,i) = u(nCells+1,i);
	u(nCells+3,i) = u(nCells+1, i);
      }

    // convert conservative values to primitive values
    u_prim = con_to_prim(u);

    // boundary extrapolation loop
    for (int i=0; i<nCells+2; i++)
      {
	for (int j=0; j<u.cols(); j++)
	  {
	    double r = slope_limiter(u(i+1, j), u(i+2, j), u(i, j));
	    double Xi = minbee(r);
	    double deltai = deltai_func(u(i+1, j), u(i+2, j), u(i, j));
	    uL(i, j) = u(i+1, j) - 0.5 * Xi * deltai;
	    uR(i, j) = u(i+1, j) + 0.5 * Xi * deltai;
	  }
      }

    uL_prim = con_to_prim(uL);
    uR_prim = con_to_prim(uR);

    uL_flux = Euler_flux_fn(uL, uL_prim);
    uR_flux = Euler_flux_fn(uR, uR_prim);

    // half-time step update
    uLhalf = uL - 0.5 * (dt / dx) * (uR_flux - uL_flux);
    uRhalf = uR - 0.5 * (dt / dx) * (uR_flux - uL_flux);
    
    uLhalf_prim = con_to_prim(uLhalf);
    uRhalf_prim = con_to_prim(uRhalf);
    
    flux = FORCE_flux(uLhalf, uRhalf, uLhalf_prim, uRhalf_prim, dt);
    //std::cout<<flux<<std::endl;
    
    // solution update loop
    for (int i = 2; i<nCells+2; i++)
      {
	for (int j=0; j<u.cols(); j++)
	  {
	    uPlus1(i, j) = u(i, j) - (dt/dx) * (flux(i-1, j) - flux(i-2, j));
	  }
      }
    u = uPlus1;
  } while (t < tStop);
}
