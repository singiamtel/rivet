// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class DM1_1978_I134061 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DM1_1978_I134061);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_npion,1,1,1);
      vector<string> energies({"483.0", "602.4", "700.0", "750.2", "773.2", "778.2", "783.2", "788.2",
          "793.2", "843.2", "893.2", "963.0", "985.0", "1003.0", "1053.0", "1096.0"});
      bool matched=false;
      for(const string& en : energies) {
        double end = std::stod(en)*MeV;
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


  RIVET_DECLARE_PLUGIN(DM1_1978_I134061);


}
