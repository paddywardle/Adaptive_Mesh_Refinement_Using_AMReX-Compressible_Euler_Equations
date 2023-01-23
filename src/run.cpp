#include "header.H"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>

void Euler1D::run()
{
  resize_matrix();
    
  initial_conds();

  u = prim_to_con(u_prim);

  uPlus1 = prim_to_con(u_prim);
  
  solvers();

  u_prim = con_to_prim(u);
}
