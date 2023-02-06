#include "SettingsData.h"

SettingsData::SettingsData()
{
  readSettings();
}

void SettingsData::readSettings()
{
  try {
    
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
    cfd_selection.lookupValue("Limiter", Limiter);
    
  }
  
  catch(...)
    {
      std::cout<<"No SettingsFile.txt present in directory."<<std::endl;
    }
}
