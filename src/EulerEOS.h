#include "common.h"

class EulerEOS
{
public:
  
  EulerEOS(double);

  Eigen::ArrayXXf prim_to_con(Eigen::ArrayXXf);

  Eigen::ArrayXXf con_to_prim(Eigen::ArrayXXf);
  
  double flux_fn_rho(double, double);

  double flux_fn_mom(double, double, double);

  double flux_fn_E(double, double, double);

  Eigen::ArrayXXf Euler_flux_fn(Eigen::ArrayXXf, Eigen::ArrayXXf);

private:

  double gamma;
    
};
