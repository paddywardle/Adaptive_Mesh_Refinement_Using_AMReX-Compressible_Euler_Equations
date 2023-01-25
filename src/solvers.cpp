#include "header.H"

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

double Euler1D::half_time_step_flux(Eigen::MatrixXd u_i, Eigen::MatrixXd u_prim_i, int col)
{
  
  double f_u;
  double rho = u_i(0);
  double energy = u_i(2);
  double velocity = u_prim_i(1);
  double pressure = u_prim_i(2);
  
  if (col == 0)
    {
      f_u = flux_fn_rho(rho, velocity);
    }
  else if (col == 1)
    {
      f_u = flux_fn_mom(rho, velocity, pressure);
    }
  else if (col == 2)
    {
      f_u = flux_fn_E(energy, velocity, pressure);
    }

  return f_u;

}

double Euler1D::calculate_timestep()
{
  
  std::vector<double> wave_speed;
  wave_speed.resize(nCells+2);

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

double Euler1D::deltai_func(double u_i, double u_iPlus1, double u_iMinus1, double w=0)
{
  // calculates cell delta value for boundary extrapolated reconstruction
  return 0.5 * (1 + w) * deltai_half(u_iMinus1, u_i); + 0.5 * (1 - w) * deltai_half(u_i, u_iPlus1);
}

double Euler1D::deltai_half(double u_i, double u_iPlus1)
{
  // calculates delta values for vertices
  return u_iPlus1 - u_i;
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
  //std::cout<<u_i<<" "<<u_iPlus1<<" "<<u_iMinus1<<std::endl;
  if (deltai_half(u_i, u_iPlus1) == 0){
    return 0;
  }
  return deltai_half(u_iMinus1, u_i) / deltai_half(u_i, u_iPlus1);
}

double Euler1D::reconstruction_XiL(double r)
{
  // calculates left slope limiter value
  return (2.0 * r) / (1 + r);
}

double Euler1D::reconstruction_XiR(double r)
{
  // calculates right slope limiter value
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

double Euler1D::lax_friedrich_flux(Eigen::MatrixXd u_i, Eigen::MatrixXd u_iPlus1, Eigen::MatrixXd u_prim_i, Eigen::MatrixXd u_prim_iPlus1, int col, double dt)
{
  double fhalf;
  
  if (col == 0)
    {
      fhalf = 0.5 * (dx / dt) * (u_i(0) - u_iPlus1(0)) + 0.5 * (flux_fn_rho(u_iPlus1(0), u_prim_iPlus1(0)) + flux_fn_rho(u_i(0), u_prim_i(0)));
    }
  else if (col == 1)
    {
      fhalf = 0.5 * (dx / dt) * (u_i(1) - u_iPlus1(1)) + 0.5 * (flux_fn_mom(u_iPlus1(0), u_prim_iPlus1(1), u_prim_iPlus1(2)) + flux_fn_mom(u_i(0), u_prim_i(1), u_prim_i(2)));
    }
  else if (col == 2)
    {
      fhalf = 0.5 * (dx / dt) * (u_i(2) - u_iPlus1(2)) + 0.5 * (flux_fn_E(u_iPlus1(2), u_prim_iPlus1(1), u_prim_iPlus1(2)) + flux_fn_E(u_i(2), u_prim_i(1), u_prim_i(2)));
    }

  return fhalf;
}

double Euler1D::richtmyer_flux(Eigen::MatrixXd u_i, Eigen::MatrixXd u_iPlus1, Eigen::MatrixXd u_prim_i, Eigen::MatrixXd u_prim_iPlus1, int col, double dt)
{
  double uhalf_rho = 0.5 * (u_i(0) + u_iPlus1(0)) - 0.5 * (dt / dx) * (flux_fn_rho(u_iPlus1(0), u_prim_iPlus1(1)) - flux_fn_rho(u_i(0), u_prim_i(1)));

  double uhalf_mom = 0.5 * (u_i(1) + u_iPlus1(1)) - 0.5 * (dt / dx) * (flux_fn_mom(u_iPlus1(0), u_prim_iPlus1(1), u_prim_iPlus1(2)) - flux_fn_mom(u_i(0), u_prim_i(1), u_prim_i(2)));

  double uhalf_E = 0.5 * (u_i(2) + u_iPlus1(2)) - 0.5 * (dt / dx) * (flux_fn_E(u_iPlus1(2), u_prim_iPlus1(1), u_prim_iPlus1(2)) - flux_fn_E(u_i(2), u_prim_i(1), u_prim_i(2)));

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

double Euler1D::FORCE_flux(Eigen::MatrixXd u_i, Eigen::MatrixXd u_iPlus1, Eigen::MatrixXd u_prim_i, Eigen::MatrixXd u_prim_iPlus1, int col, double dt)
{
  double fhalf = 0.5 * (lax_friedrich_flux(u_i, u_iPlus1, u_prim_i, u_prim_iPlus1, col, dt) + richtmyer_flux(u_i, u_iPlus1, u_prim_i, u_prim_iPlus1, col, dt));

  return fhalf;
}

void Euler1D::solvers()
{
  // sets initial conditions into instance variables
  initial_conds();

  // set current time to simulation start time
  double t = tStart;

  Eigen::MatrixXd flux(nCells+1, 3);
  Eigen::MatrixXd uL(nCells, 3);
  Eigen::MatrixXd uR(nCells, 3);
  Eigen::MatrixXd uL_prim(nCells, 3);
  Eigen::MatrixXd uR_prim(nCells, 3);
  Eigen::MatrixXd uLhalf(nCells, 3);
  Eigen::MatrixXd uRhalf(nCells, 3);
  Eigen::MatrixXd uLhalf_prim(nCells, 3);
  Eigen::MatrixXd uRhalf_prim(nCells, 3);
  //std::cout<<u<<std::endl;
  
  do {
    // calculate timestep based on initial data
    double dt = calculate_timestep();

    // convert conservative values to primitive values
    u_prim = con_to_prim(u);

    // add timestep to current time
    t += dt;

    // add transmissive boundary
    for (int i=0; i<u.cols(); i++)
      {
	u(0,i) = u(1, i);
	u(nCells+1, i) = u(nCells, i);
      }

    //std::cout<<u<<std::endl;

    // boundary extrapolation loop
    for (int i=1; i<nCells; i++)
      {
	for (int j=0; j<u.cols(); j++)
	  {
	    double r = slope_limiter(u(i, j), u(i+1, j), u(i-1, j));
	    double XiL = reconstruction_XiL(r);
	    double XiR = reconstruction_XiR(r);
	    double Xi = minibee(r, XiR);
	    double deltai = deltai_func(u(i, j), u(i+1, j), u(i-1, j));
	    uL(i, j) = reconstruction_uL(u(i, j), Xi, deltai);
	    uR(i, j) = reconstruction_uR(u(i, j), Xi, deltai);
	    //std::cout<<u(i, j)<<" "<<u(i+1, j)<<" "<<u(i-1, j)<<std::endl;
	    std::cout<<r<<" "<<XiL<<" "<<XiR<<" "<<Xi<<" "<<deltai<<std::endl;
	  }
      }
    
    uL_prim = con_to_prim(uL);
    uR_prim = con_to_prim(uR);

    // half-time step update
    for (int i=0; i<uL.rows(); i++)
      {
	for (int j=0; j<uL.cols(); j++)
	  {
	    uLhalf(i, j) = uL(i, j) - 0.5 * (dt / dx) * (half_time_step_flux(uR.row(i), uR_prim.row(i), j) - half_time_step_flux(uL.row(i), uL_prim.row(i), j));
	    uRhalf(i, j) = uR(i, j) - 0.5 * (dt / dx) * (half_time_step_flux(uR.row(i), uR_prim.row(i), j) - half_time_step_flux(uL.row(i), uL_prim.row(i), j));
	  }
      }

    uLhalf_prim = con_to_prim(uLhalf);
    uRhalf_prim = con_to_prim(uRhalf);

    // flux calculation
    for (int i=0; i<nCells; i++)
      {
	for (int j=0; j<flux.cols(); j++)
	  {
	    flux(i, j) = FORCE_flux(uLhalf.row(i), uRhalf.row(i), uLhalf_prim.row(i), uRhalf_prim.row(i), j, dt);
	  }
      }
    
    // solution update loop
    for (int i = 1; i<nCells+1; i++)
      {
	for (int j=0; j<u.cols(); j++)
	  {
	    uPlus1(i, j) = u(i, j) - (dt/dx) * (flux(i, j) - flux(i-1, j));
	  }
      }
    u = uPlus1;
    break;
  } while (t < tStop);
}
