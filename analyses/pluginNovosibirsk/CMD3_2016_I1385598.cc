// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief e+ e- > p pbar
  class CMD3_2016_I1385598 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD3_2016_I1385598);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_nproton, 1,1,6);
      for (const string& en : _nproton.binning().edges<0>()) {
        double end = en=="1900 (2012)"s ? 1.9 : std::stod(en)*MeV;
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
      if (fs.particles().size() != 2) vetoEvent;
      for (const Particle& p : fs.particles()) {
        if (p.abspid() != PID::PROTON) vetoEvent;
      }
      _nproton->fill(_ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nproton, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nproton;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD3_2016_I1385598);


}
