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
  double nCells;
  double dx = (x1-x0)/nCells;

  Eigen::MatrixXd u;

  Eigen::MatrixXd uPlus1;

  Eigen::MatrixXd u_prim;

  // private member functions

  void resize_matrix();

  void initial_conds();

  Eigen::MatrixXd prim_to_con(Eigen::MatrixXd);

  Eigen::MatrixXd con_to_prim(Eigen::MatrixXd);

  double deltai_func(double, double, double, double);

  double deltai_half(double, double);
  
  double reconstruction_uL(double, double, double);

  double reconstruction_uR(double, double, double);

  double slope_limiter(double, double, double);

  double reconstruction_XiL(double);

  double reconstruction_XiR(double);

  double minibee(double, double);
  
  double calculate_timestep();
  
  double flux_fn_rho(double, double);

  double flux_fn_mom(double, double, double);

  double flux_fn_E(double, double, double);

  double lax_friedrich_flux(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, int, double);

  double richtmyer_flux(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, int, double);

  double half_time_step_flux(Eigen::MatrixXd, Eigen::MatrixXd, int);

  double FORCE_flux(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, int, double);

  void solvers();
    
};
