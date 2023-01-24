#include <vector>
#include <array>
#include <iostream>
#include "header.H"

int main()
{
  double x0 = 0.0;
  double x1 = 1.0;
  double tStart = 0.0;
  double tStop = 0.25;
  double gamma = 1.4;
  double C = 0.8;
  int nCells = 100;

  std::cin>>x0>>x1>>tStart>>tStop>>gamma>>C>>nCells;

  Euler1D E(nCells, tStart, tStop, x0, x1, gamma, C);

  E.run();

  E.outputFile("../data/results.dat");
}
