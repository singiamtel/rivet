// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Cross section for $e^+e^-\to$ proton and antiproton between 2. and 3.8 GeV using ISR
  class BESIII_2019_I1718337 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1718337);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_nproton, 1,1,1);
      for (const string& en : _nproton.binning().edges<0>()) {
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
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
        if(abs(p.pid())!=PID::PROTON) vetoEvent;
      }
      _nproton->fill(_ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nproton, crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nproton;
    string _ecms;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1718337);

}
