// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Measurement of the hadronic cross section for energies between 3.08 and 3.116 GeV
  class BES_1995_I39870 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BES_1995_I39870);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons, 1, 1, 1);
      for (const string& en : _c_hadrons.binning().edges<0>()) {
        double end = std::stod(en)*GeV;
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
      // mu+mu- + photons
      if (nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) vetoEvent;
      // everything else
      else _c_hadrons->fill(_ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      scale(_c_hadrons,fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons;
    string _ecms;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BES_1995_I39870);

}
