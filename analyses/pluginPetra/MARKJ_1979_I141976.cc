// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MARKJ_1979_I141976 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKJ_1979_I141976);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons, "sigma_hadrons", refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_muons, "sigma_muons", refData<YODA::BinnedEstimate<string>>(1,1,1));
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
      if (nCount[-13]==1 and nCount[13]==1 &&  ntotal==2+nCount[22]) {
        _c_muons->fill(Ecm);
      }
      // everything else
      else {
        _c_hadrons->fill(Ecm);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/ sumOfWeights() /picobarn;
      scale({_c_hadrons, _c_muons}, fact);
      BinnedEstimatePtr<string> mult;
      book(mult, 1, 1, 1);
      divide(_c_hadrons, _c_muons, mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons, _c_muons;
    const string Ecm = "31.57";
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MARKJ_1979_I141976);


}
