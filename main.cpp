#include "simulation_x.hpp"
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cout << "Wrong number of arguments supplied" << std::endl;
    return 0;
  }
  std::ifstream in_file;
  in_file.open(std::string(argv[1]));
  std::ofstream out_file;
  out_file.open(std::string(argv[2]));
  XSimulation simul;
  /*simul.Initial(in_file);
  simul.Dump(out_file);
  simul.FillENk();
  simul.FillVB();
  simul.FillForceVB();
  simul.DumpForce(out_file);
  simul.DumpEnergy(out_file);
  simul.UpdateMNHC_VV3();
  simul.UpdateMNHC_VV3();
  //simul.Dump(out_file);
  for (int i = 0; i < 1000; ++i)
    simul.UpdateMNHC_VV3();
  simul.Dump(out_file);*/
  simul.Initial(in_file);
  simul.even = true;
  //simul.SetVi(1.0);
  XNum angle = simul.GetVi()*10*(M_PI/2.0/20.0);
  simul.SetVi(1.0*std::cos(angle));
  XSimulation simul2;
  in_file.close();
  in_file.open(std::string(argv[1]));
  simul2.Initial(in_file);
  simul2.even = false;
  simul2.SetVi(1.0*std::sin(angle));
  simul.vi2 = simul2.GetVi();
  simul2.vi2 = simul.GetVi();
  int i;
  clock_t t;
  t = clock();
  for (i = 0; i < simul.skip; ++i) {
    simul.PeriodBoundary();
    simul.UpdateMNHC_VV3();
    if (!simul.ok) {
      simul.Dump(out_file);
      simul.DumpForce(out_file);
      simul.DumpEnergy(out_file);
      simul.NormalizeStat();
      simul.DumpStat(out_file);
      return 0;
    }
    simul2.PeriodBoundary();
    simul2.UpdateMNHC_VV3();
    if (!simul2.ok) {
      simul2.Dump(out_file);
      simul2.DumpForce(out_file);
      simul2.DumpEnergy(out_file);
      simul2.NormalizeStat();
      simul2.DumpStat(out_file);
      return 0;
    }
  }
  for (i = 0; i < simul.step; ++i) {
    simul.PeriodBoundary();
    simul.UpdateMNHC_VV3();
    if (!simul.ok) {
      simul.Dump(out_file);
      simul.DumpForce(out_file);
      simul.DumpEnergy(out_file);
      simul.NormalizeStat();
      simul.DumpStat(out_file);
      return 0;
    }
    simul2.PeriodBoundary();
    simul2.UpdateMNHC_VV3();
    if (!simul2.ok) {
      simul2.Dump(out_file);
      simul2.DumpForce(out_file);
      simul2.DumpEnergy(out_file);
      simul2.NormalizeStat();
      simul2.DumpStat(out_file);
      return 0;
    }
    simul.MakeStat();
    simul2.MakeStat();
    if (i % 100000 == 0)
      out_file << "i=" << i << std::endl;
  }
  //simul.Dump(out_file);
  simul.NormalizeStat();
  simul.DumpStat(out_file);
  simul2.NormalizeStat();
  simul2.DumpStat(out_file);
  t = clock() - t;
  out_file << (int)(((double)1000 * t) / CLOCKS_PER_SEC) << std::endl;
  return 0;
}
