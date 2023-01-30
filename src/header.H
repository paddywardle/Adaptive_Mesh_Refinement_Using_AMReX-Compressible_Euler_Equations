#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <libconfig.h++>
#include <eigen3/Eigen/Dense>

class Euler1D
{
public:
  Euler1D(int, double, double, double, double, double, double);

  void outputFile(std::string);

  void run();

private:
  
  // member data
  
  double tStart;
  double tStop;
  double x0;
  double x1;
  double gamma;
  double C;
  int nCells;
  double dx = (x1-x0)/nCells;

  Eigen::ArrayXXf u;

  Eigen::ArrayXXf uPlus1;

  Eigen::ArrayXXf u_prim;

  // private member functions

  void resize_matrix();

  void initial_conds();

  Eigen::ArrayXXf prim_to_con(Eigen::ArrayXXf);

  Eigen::ArrayXXf con_to_prim(Eigen::ArrayXXf);

  double deltai_func(double, double, double, double);
  
  double reconstruction_uL(double, double, double);

  double reconstruction_uR(double, double, double);

  double slope_limiter(double, double, double);

  double reconstruction_XiR(double);

  double minbee(double);
  
  double calculate_timestep();
  
  double flux_fn_rho(double, double);

  double flux_fn_mom(double, double, double);

  double flux_fn_E(double, double, double);

  Eigen::ArrayXXf Euler_flux_fn(Eigen::ArrayXXf, Eigen::ArrayXXf);

  double lax_friedrich_flux(double, double, double, double, double);

  double richtmyer_flux(double, double, double, double, double);

  Eigen::ArrayXXf FORCE_flux(Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf, double);

  void solvers();
    
};
