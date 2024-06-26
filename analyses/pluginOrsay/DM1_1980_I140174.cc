// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > pi+ pi- pi0
  class DM1_1980_I140174 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DM1_1980_I140174);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_n3pi, 1, 1, 1);
      for (const string& en : _n3pi.binning().edges<0>()) {
        double end = std::stod(en)*MeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      if(ntotal==3 && nCount[211] == 1 && nCount[-211] == 1 && nCount[111] == 1)
	_n3pi->fill(_ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_n3pi, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _n3pi;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DM1_1980_I140174);


}
