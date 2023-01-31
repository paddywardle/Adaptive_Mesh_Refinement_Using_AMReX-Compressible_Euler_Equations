#include "EulerSystem.h"

int main()
{
  
  EulerSystem E;
  
  E.run();

  E.outputFile("../data/results.dat");
}
