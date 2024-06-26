// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief e+e- > 2p 2pbar
  class BESIII_2020_I1837725 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1837725);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_nproton, 1, 1, 1);
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
      if(fs.particles().size()!=4) vetoEvent;
      for (const Particle& p :  fs.particles()) {
	if(abs(p.pid())!=PID::PROTON) vetoEvent;
      }
      _nproton->fill(_ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nproton, crossSection()/ sumOfWeights() /femtobarn);
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _nproton;
    string _ecms;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1837725);

}
