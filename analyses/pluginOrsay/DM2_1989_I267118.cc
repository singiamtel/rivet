// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class DM2_1989_I267118 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DM2_1989_I267118);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_npion,1,1,1);
      vector<string> energies({"1.35", "1.395", "1.425", "1.475", "1.525", "1.575", "1.625", "1.675", "1.725",
          "1.775", "1.825", "1.875", "1.925", "1.975", "2.025", "2.075", "2.125"});
      bool matched=false;
      for(const string& en : energies) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          ecms = en;
          matched=true;
          break;
        }
      }
      if(!matched) MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
	if(abs(p.pid())!=PID::PIPLUS) vetoEvent;
      }
      _npion->fill(ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_npion,crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _npion;
    string ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DM2_1989_I267118);


}
