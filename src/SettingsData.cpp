#include "SettingsData.h"

SettingsData::SettingsData()
{
  readSettings();
}

void SettingsData::readSettings()
{
  try {

    std::string limiter;
    
    libconfig::Config cfg;
    cfg.readFile("SettingsFile.txt");
    const auto& cfd_selection = cfg.lookup("CFD");
    const auto& spatial_domain = cfd_selection.lookup("SpatialDomain");
    const auto& temporal_domain = cfd_selection.lookup("TemporalDomain");
    const auto& EOS = cfd_selection.lookup("EOS");

    spatial_domain.lookupValue("x0", x0);
    spatial_domain.lookupValue("x1", x1);
    spatial_domain.lookupValue("nCells", nCells);
    temporal_domain.lookupValue("tStart", tStart);
    temporal_domain.lookupValue("tStop", tStop);
    EOS.lookupValue("gamma", gamma);
    cfd_selection.lookupValue("C", C);
    cfd_selection.lookupValue("Limiter", limiter);

    if (limiter == "Superbee")
      {
	Lim = Limiters::Superbee;
      }
    else if (limiter == "Van_Leer")
      {
	Lim = Limiters::Van_Leer;
      }
    else if (limiter == "Van_Albada")
      {
	Lim = Limiters::Van_Albada;
      }
    else if (limiter == "Minbee")
      {
	Lim = Limiters::Minbee;
      }
  }
  
  catch(...)
    {
      std::cout<<"No SettingsFile.txt present in directory."<<std::endl;
    }
}
