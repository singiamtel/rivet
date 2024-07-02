// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief 
  class BESIII_2020_I1808875 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1808875);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Projections
      declare(FinalState(), "FS");
      book(_sigma_muons, 1,1,1);
      for (const string& en : _sigma_muons.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
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
      // mu+mu- + photons
      if (nCount[-13]==1 && nCount[13]==1 && ntotal==2+nCount[22]) {
        _sigma_muons->fill(_ecms);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma_muons, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma_muons;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1808875);

}
