#include <libconfig.h++>
#include <iostream>

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
  
};
