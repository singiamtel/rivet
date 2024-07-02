// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Charged multiplicity
  class MARKI_1982_I169326 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKI_1982_I169326);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(ChargedFinalState(), "FS");

      // Book histograms
      book(_nHadrons, 6, 1, 1);
      for (const string& en : _nHadrons.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      _nHadrons->fill(_ecms,fs.particles().size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {}
    /// @}


    /// @name Histograms
    /// @{
    BinnedProfilePtr<string> _nHadrons;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(MARKI_1982_I169326);


}
