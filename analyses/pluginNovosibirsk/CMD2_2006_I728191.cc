// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class CMD2_2006_I728191 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD2_2006_I728191);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_npion, 3,1,1);
      vector<int> enint({370, 390, 410, 430, 450, 470, 480, 500, 510, 520});
      icms=-1;
      for(const int& en : enint) {
        double end = double(en)*MeV;
        if(isCompatibleWithSqrtS(end)) {
          icms = en;
          break;
        }
      }
      if(icms<0) MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
	if(abs(p.pid())!=PID::PIPLUS) vetoEvent;
      }
      _npion->fill(icms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_npion,crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _npion;
    int icms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD2_2006_I728191);


}
