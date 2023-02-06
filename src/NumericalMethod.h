#include "EulerEOS.h"
#include "Limiters.h"

class NumericalMethod
{
public:

  NumericalMethod(double, int);

  double reconstruction_uL(double, double, double, Limiters);

  double reconstruction_uR(double, double, double, Limiters);

  Eigen::ArrayXXf FORCE_flux(Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf, double, double);

  Eigen::ArrayXXf HLLC_flux(Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf, double, double);

  Eigen::ArrayXXf uL_half_update(Eigen::ArrayXXf, Eigen::ArrayXXf, double, double);

  Eigen::ArrayXXf uR_half_update(Eigen::ArrayXXf, Eigen::ArrayXXf, double, double);

private:

  double gamma;
  int nCells;
  EulerEOS eos;

  double deltai_func(double, double, double, double);

  double slope_limiter(double, double, double);

  double minbee(double);

  double superbee(double);

  double van_leer(double);

  double van_albada(double);

  double limiterXi(double, Limiters);

  double lax_friedrich_flux(double, double, double, double, double, double);

  double richtmyer_flux(double, double, double, double, double, double);

  Eigen::ArrayXXf wavespeed(Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf);

  Eigen::ArrayXXf uHLLC(Eigen::ArrayXXf, Eigen::ArrayXXf, double, double);

  Eigen::ArrayXXf fHLLC(Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf, double, double, double);

  
  
};
