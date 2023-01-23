#include <vector>
#include <string>
#include <array>

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

  std::vector<std::array<double, 3>> u;

  std::vector<std::array<double, 3>> uPlus1;

  std::vector<std::array<double, 3>> u_prim;

  // private member functions

  void resize_matrix();

  void initial_conds();

  std::vector<std::array<double, 3>> prim_to_con(std::vector<std::array<double, 3>>);

  std::vector<std::array<double, 3>> con_to_prim(std::vector<std::array<double, 3>>);

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

  double lax_friedrich_flux(std::array<double, 3>, std::array<double, 3>, std::array<double, 3>, std::array<double, 3>, int, double);

  double richtmyer_flux(std::array<double, 3>, std::array<double, 3>, std::array<double, 3>, std::array<double, 3>, int, double);

  double half_time_step_flux(std::array<double, 3>, std::array<double, 3>, int);

  double FORCE_flux(std::array<double, 3>, std::array<double, 3>, std::array<double, 3>, std::array<double, 3>, int, double);

  void solvers();
    
};
