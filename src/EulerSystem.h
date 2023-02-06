#include "SettingsData.h"
#include "NumericalMethod.h"
#include "common.h"

class EulerSystem: public SettingsData
{
public:
  
  EulerSystem();

  void outputFile(std::string);

  void run();

private:
  
  double dx = (x1-x0)/nCells;

  Eigen::ArrayXXf u;
  Eigen::ArrayXXf uPlus1;
  Eigen::ArrayXXf u_prim;
  Limiters Lim;

  // private member functions

  void resize_matrix();

  void initial_conds();
  
  double calculate_timestep();

  void Solver1D();
};
