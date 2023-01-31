#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <libconfig.h++>
#include <eigen3/Eigen/Dense>
#include "EulerEOS.h"

class NumericalMethod
{
public:

  NumericalMethod(double, int);

  double reconstruction_uL(double, double, double);

  double reconstruction_uR(double, double, double);

  Eigen::ArrayXXf FORCE_flux(Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf, double, double);

  Eigen::ArrayXXf uL_half_update(Eigen::ArrayXXf, Eigen::ArrayXXf, double, double);

  Eigen::ArrayXXf uR_half_update(Eigen::ArrayXXf, Eigen::ArrayXXf, double, double);

private:

  double gamma;
  int nCells;
  EulerEOS eos;

  double deltai_func(double, double, double, double);

  double slope_limiter(double, double, double);

  double minbee(double);

  double lax_friedrich_flux(double, double, double, double, double, double);

  double richtmyer_flux(double, double, double, double, double, double);
  
};
