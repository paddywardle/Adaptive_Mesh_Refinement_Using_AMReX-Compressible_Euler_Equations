#include "header.H"

Eigen::MatrixXd Euler1D::Euler_flux_fn(Eigen::MatrixXd f, Eigen::MatrixXd f_prim)
{
  
  Eigen::MatrixXd flux_fn(f.rows(), f.cols());
  
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

double Euler1D::deltai_func(double u_i, double u_iPlus1, double u_iMinus1, double w=0)
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
  //std::cout<<u_i<<" "<<u_iPlus1<<" "<<u_iMinus1<<std::endl;
  if ((u_iPlus1 - u_i) == 0.0){
    return 0.0;
  }
  return (u_i - u_iMinus1) / (u_iPlus1 - u_i);
}

double Euler1D::reconstruction_XiR(double r)
{
  // calculates right slope limiter value
  return 2.0 / (1.0 + r);
}

double Euler1D::minbee(double r, double XiR)
{
  if (r <= 0.0)
    {
      return 0.0;
    }
  else if (r > 1.0)
    {
      return std::min(1.0, XiR);
    }
  else
    {
      return r;
    }
}

Eigen::MatrixXd Euler1D::lax_friedrich_flux(Eigen::MatrixXd u_i, Eigen::MatrixXd u_iPlus1,  Eigen::MatrixXd u_prim_i, Eigen::MatrixXd u_prim_iPlus1, double dt)
{

  Eigen::MatrixXd fhalf(u_i.rows(), u_i.cols());

  fhalf = 0.5 * (dx / dt) * (u_i - u_iPlus1) + 0.5 * (Euler_flux_fn(u_iPlus1, u_prim_iPlus1) - Euler_flux_fn(u_i, u_prim_i));
  
  return fhalf;
}

Eigen::MatrixXd Euler1D::richtmyer_flux(Eigen::MatrixXd u_i, Eigen::MatrixXd u_iPlus1, Eigen::MatrixXd u_prim_i, Eigen::MatrixXd u_prim_iPlus1,  double dt)
{

  Eigen::MatrixXd uhalf(u_i.rows(), u_i.cols());
  Eigen::MatrixXd fhalf(u_i.rows(), u_i.cols());

  uhalf = 0.5 * (u_i - u_iPlus1) - 0.5 * (dt / dx) * (Euler_flux_fn(u_iPlus1, u_prim_iPlus1) - Euler_flux_fn(u_i, u_prim_i));

  fhalf = Euler_flux_fn(uhalf, con_to_prim(uhalf));

  return fhalf;
}

Eigen::MatrixXd Euler1D::FORCE_flux(Eigen::MatrixXd u_i, Eigen::MatrixXd u_iPlus1, Eigen::MatrixXd u_prim_i, Eigen::MatrixXd u_prim_iPlus1, double dt)
{
  Eigen::MatrixXd fhalf(u_i.rows(), u_i.cols());

  fhalf = lax_friedrich_flux(u_i, u_iPlus1, u_prim_i, u_prim_iPlus1, dt);
  std::cout<<fhalf.rows()<<" "<<fhalf.cols()<<std::endl;

  fhalf = richtmyer_flux(u_i, u_iPlus1, u_prim_i, u_prim_iPlus1, dt);
  std::cout<<fhalf.rows()<<" "<<fhalf.cols()<<std::endl;

  fhalf = 0.5 * (lax_friedrich_flux(u_i, u_iPlus1, u_prim_i, u_prim_iPlus1, dt) + richtmyer_flux(u_i, u_iPlus1, u_prim_i, u_prim_iPlus1, dt));

  return fhalf;
}

void Euler1D::solvers()
{
  // sets initial conditions into instance variables
  initial_conds();

  // set current time to simulation start time
  double t = tStart;

  Eigen::MatrixXd flux(nCells+1, 3);
  Eigen::MatrixXd uL(nCells+2, 3);
  Eigen::MatrixXd uR(nCells+2, 3);
  Eigen::MatrixXd uL_prim(nCells+2, 3);
  Eigen::MatrixXd uR_prim(nCells+2, 3);
  Eigen::MatrixXd uLhalf(nCells+2, 3);
  Eigen::MatrixXd uRhalf(nCells+2, 3);
  Eigen::MatrixXd uLhalf_prim(nCells+2, 3);
  Eigen::MatrixXd uRhalf_prim(nCells+2, 3);

  Eigen::MatrixXd uL_flux(nCells+2, 3);
  Eigen::MatrixXd uR_flux(nCells+2, 3);
  
  do {
    
    // calculate timestep based on initial data
    double dt = calculate_timestep();
    std::cout<<t<<" "<<dt<<std::endl;
    
    // convert conservative values to primitive values
    u_prim = con_to_prim(u);

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
    
    // boundary extrapolation loop
    for (int i=0; i<nCells+2; i++)
      {
	for (int j=0; j<u.cols(); j++)
	  {
	    double r = slope_limiter(u(i+1, j), u(i+2, j), u(i, j));
	    double XiR = reconstruction_XiR(r);
	    double Xi = minbee(r, XiR);
	    double deltai = deltai_func(u(i+1, j), u(i+2, j), u(i, j));
	    uL(i, j) = reconstruction_uL(u(i+1, j), Xi, deltai);
	    uR(i, j) = reconstruction_uR(u(i+1, j), Xi, deltai);
	    //std::cout<<r<<" "<<XiR<<" "<<Xi<<" "<<deltai<<std::endl;
	  }
      }
    
    uL_prim = con_to_prim(uL);
    uR_prim = con_to_prim(uR);

    uL_flux = Euler_flux_fn(uL, uL_prim);
    uR_flux = Euler_flux_fn(uR, uR_prim);
    std::cout<<uL_prim<<std::endl;

    // half-time step update
    uLhalf = uL - 0.5 * (dt / dx) * (uR_flux - uL_flux);
    uRhalf = uR - 0.5 * (dt / dx) * (uR_flux - uL_flux);

    uLhalf_prim = con_to_prim(uLhalf);
    uRhalf_prim = con_to_prim(uRhalf);

    std::cout<<flux.rows()<<" "<<flux.cols()<<std::endl;
    flux = FORCE_flux(uLhalf, uRhalf, uLhalf, uRhalf_prim, dt);
    std::cout<<flux.rows()<<" "<<flux.cols()<<std::endl;
    
    // solution update loop
    for (int i = 2; i<nCells+2; i++)
      {
	for (int j=0; j<u.cols(); j++)
	  {
	    uPlus1(i, j) = u(i, j) - (dt/dx) * (flux(i-1, j) - flux(i-2, j));
	  }
      }
    u = uPlus1;
    break;
  } while (t < tStop);
}
