// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class JADE_1979_I142874 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(JADE_1979_I142874);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(ChargedFinalState(), "FS");
      if     ( isCompatibleWithSqrtS(22.0*GeV) ) _sqs = "22.0";
      else if( isCompatibleWithSqrtS(27.7*GeV) ) _sqs = "27.7";
      else if( isCompatibleWithSqrtS(30.0*GeV) ) _sqs = "30.0";
      else if( isCompatibleWithSqrtS(31.6*GeV) ) _sqs = "31.6";
      else
        MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
                    << " doesn't match any available analysis energy .");
      // Book histograms
      book(_mult, 2, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      _mult->fill(_sqs,fs.particles().size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
    }
    /// @}


    /// @name Histograms
    /// @{
    BinnedProfilePtr<string> _mult;
    string _sqs;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(JADE_1979_I142874);


}
