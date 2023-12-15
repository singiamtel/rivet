// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class VENUS_1999_I500179 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(VENUS_1999_I500179);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons, 1, 1, 1);
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
        vetoEvent; // mu+mu- + photons
      }
      else {
        _c_hadrons->fill(Ecm); // everything else
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/sumOfWeights()/picobarn;
      scale(_c_hadrons, fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons;
    const string Ecm = "57.77";
    /// @}


  };


  RIVET_DECLARE_PLUGIN(VENUS_1999_I500179);


}
