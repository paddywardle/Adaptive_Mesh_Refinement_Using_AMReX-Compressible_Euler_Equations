#include "header.H"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>

Euler1D::Euler1D(int nCells, double tStart, double tStop, double x0, double x1, double gamma, double C)
  :nCells(nCells), tStart(tStart), tStop(tStop), x0(x0), x1(x1), gamma(gamma), C(C){}

void Euler1D::resize_matrix()
{
  u.resize(nCells+2, 3);
  uPlus1.resize(nCells+2, 3);
  u_prim.resize(nCells+2, 3);
}

