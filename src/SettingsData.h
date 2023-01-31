#include <libconfig.h++>
#include <iostream>
#include <string>
#include "Limiters.h"

class SettingsData
{
 public:
  
  SettingsData();

  void readSettings();
  
 protected:

  double x0;
  double x1;
  double tStart;
  double tStop;
  double gamma;
  double C;
  int nCells;
  Limiters Lim;
  
};
