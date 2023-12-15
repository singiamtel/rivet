// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class NMD_1974_I745 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(NMD_1974_I745);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons, 1, 1, 1);
      book(_c_muons, "sigma_muons", refData<YODA::BinnedEstimate<int>>(1,1,1));
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
      if (nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) {
        _c_muons->fill(4); // mu+mu- + photons
      }
      else {
        _c_hadrons->fill(4); // everything else
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/sumOfWeights()/nanobarn;
      scale( {_c_hadrons, _c_muons}, fact);
      BinnedEstimatePtr<int> mult;
      book(mult, 1, 1, 2);
      divide(_c_hadrons, _c_muons, mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _c_hadrons, _c_muons;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(NMD_1974_I745);


}
