// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- > pi+pi-
  class KLOE_2005_I655225 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(KLOE_2005_I655225);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // Book histograms
      book(_npion, 2, 1, 1);
      for (const string& en : _npion.binning().edges<0>()) {
        double end = sqrt(std::stod(en))*GeV;
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
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
	if(abs(p.pid())!=PID::PIPLUS) vetoEvent;
      }
      _npion->fill(_ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_npion,crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string>  _npion;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(KLOE_2005_I655225);


}
