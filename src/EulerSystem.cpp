#include "EulerSystem.h"

EulerSystem::EulerSystem()
  :SettingsData(){};

double EulerSystem::calculate_timestep()
{
  EulerEOS eos(gamma);
  
  std::vector<double> wave_speed;
  wave_speed.resize(nCells+4);

  u_prim = eos.con_to_prim(u);
  
  for (int i=0; i<u.rows(); i++)
    {
      double cs = sqrt((gamma*u_prim(i, 2))/u_prim(i, 0));
      wave_speed[i] = abs(u_prim(i, 1)) + cs;
    }

  double a_max = *std::max_element(wave_speed.begin(), wave_speed.end());

  double dt = C * (dx / a_max);

  return dt;
}

void EulerSystem::initial_conds()
{
  for (int i=0; i<u.rows(); i++)
    {
      double x = x0 + i * dx;
      
      if (x <= 0.5)
	{
	  // density
	  u_prim(i, 0) = 1.0;
	  // velocity
	  u_prim(i, 1) = 0.0;
	  // pressure
	  u_prim(i, 2) = 1.0;
	}
      else
	{
	  // density
	  u_prim(i, 0) = 0.125;
	  // velocity
	  u_prim(i, 1) = 0.0;
	  // pressure
	  u_prim(i, 2) = 0.1;
	}
    }
}

void EulerSystem::resize_matrix()
{
  u.resize(nCells+4, 3);
  uPlus1.resize(nCells+4, 3);
  u_prim.resize(nCells+4, 3);
}

void EulerSystem::outputFile(std::string outputName)
{
  std::ofstream output(outputName);

  for (int i=1; i<u_prim.rows()-2; i++)
    {
      double x = x0 + (i-1) * dx;
      output<<x<<" ";
      for (int j=0; j<u_prim.cols(); j++)
	{
	  output<<u_prim(i, j)<<" ";
	}
      output<<std::endl;
    }
}

void EulerSystem::Solver1D()
{
  // sets initial conditions into instance variables
  initial_conds();

  // set current time to simulation start time
  double t = tStart;

  // Instance of helper classes
  EulerEOS eos(gamma);
  NumericalMethod NumM(gamma, nCells);

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
    u_prim = eos.con_to_prim(u);

    // boundary extrapolation loop
    for (int i=0; i<nCells+2; i++)
      {
	for (int j=0; j<u.cols(); j++)
	  {
	    uL(i, j) = NumM.reconstruction_uL(u(i+1, j), u(i+2, j), u(i, j), Lim);
	    uR(i, j) = NumM.reconstruction_uR(u(i+1, j), u(i+2, j), u(i, j), Lim);
	  }
      }

    // half-time step update
    uLhalf = NumM.uL_half_update(uL, uR, dt, dx);
    uRhalf = NumM.uR_half_update(uL, uR, dt, dx);
    
    uLhalf_prim = eos.con_to_prim(uLhalf);
    uRhalf_prim = eos.con_to_prim(uRhalf);
    
    flux = NumM.FORCE_flux(uLhalf, uRhalf, uLhalf_prim, uRhalf_prim, dt, dx);
    
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

void EulerSystem::run()
{

  EulerEOS eos(gamma);
  
  resize_matrix();
    
  initial_conds();
  
  u = eos.prim_to_con(u_prim);
  
  uPlus1 = eos.prim_to_con(u_prim);
  
  Solver1D();

  u_prim = eos.con_to_prim(u);
}
