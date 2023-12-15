// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PLUTO_1979_I142517 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PLUTO_1979_I142517);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_c_hadrons, "sigma_hadrons", refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_muons, "sigma_muons", refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_mult, 1, 1, 1);
      for (const string label : {"22.0", "27.6", "30.0", "31.6"}) {
        const double E = std::stod(label);
        if (isCompatibleWithSqrtS(E*GeV))  Ecm = label;
      }
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
        _c_muons->fill(Ecm); // mu+mu- + photons
      }
      else {
        _c_hadrons->fill(Ecm); // everything else
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/ sumOfWeights() /nanobarn;
      scale({_c_hadrons, _c_muons}, fact);
      divide(_c_hadrons, _c_muons, _mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons, _c_muons;
    BinnedEstimatePtr<string> _mult;
    string Ecm = "";
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PLUTO_1979_I142517);


}
