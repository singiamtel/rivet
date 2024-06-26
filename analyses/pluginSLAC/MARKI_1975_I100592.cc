// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief cross section and charged multiplicity
  class MARKI_1975_I100592 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKI_1975_I100592);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(ChargedFinalState(), "FS");

      // Book histograms
      book(_nEvent  , 1, 1, 1);
      book(_nHadrons, 2, 1, 1);
      for (const string& en : _nEvent.binning().edges<0>()) {
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
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      _nEvent->fill(_ecms);
      _nHadrons->fill(_ecms,fs.particles().size());
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nEvent,crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nEvent;
    BinnedProfilePtr<string> _nHadrons;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(MARKI_1975_I100592);


}
