#include "header.H"

int main()
{
  double x0;
  double x1;
  double tStart;
  double tStop;
  double gamma;
  double C;
  int nCells;

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

  std::cout<<x0<<" "<<x1<<" "<<nCells<<" "<<tStart<<" "<<tStop<<" "<<gamma<<" "<<C<<std::endl;

  Euler1D E(nCells, tStart, tStop, x0, x1, gamma, C);

  E.run();

  E.outputFile("../data/results.dat");
}
