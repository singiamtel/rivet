// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- -> pi+pi-2pi0 and 2pi0gamma near the phi
  class KLOE_2008_I791841 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(KLOE_2008_I791841);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_n4pi, 1,1,1);
      book(_n2pigamma, 2,1,1);
      for (const string& en : _n4pi.binning().edges<0>()) {
        const double end = std::stod(en)*MeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
      }
      if (nCount[111]==2) {
        if( nCount[211] == 1 && nCount[-211] == 1 ) {
          _n4pi->fill(_ecms);
        }
        else if( nCount[22] == 1) {
          _n2pigamma->fill(_ecms);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      scale(_n4pi,fact);
      scale(_n2pigamma,fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _n4pi,_n2pigamma;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(KLOE_2008_I791841);


}
