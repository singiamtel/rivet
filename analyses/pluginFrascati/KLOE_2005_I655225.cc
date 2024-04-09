// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
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
      vector<string> energies({"0.35", "0.36", "0.37", "0.38", "0.39", "0.4", "0.41", "0.42", "0.43",
          "0.44", "0.45", "0.46", "0.47", "0.48", "0.49", "0.5", "0.51", "0.52", "0.53", "0.54",
          "0.55", "0.56", "0.57", "0.58", "0.59", "0.6", "0.61", "0.62", "0.63", "0.64", "0.65",
          "0.66", "0.67", "0.68", "0.69", "0.7", "0.71", "0.72", "0.73", "0.74", "0.75", "0.76",
          "0.77", "0.78", "0.79", "0.8", "0.81", "0.82", "0.83", "0.84", "0.85", "0.86", "0.87",
          "0.88", "0.89", "0.9", "0.91", "0.92", "0.93", "0.94"});
      bool matched=false;
      for(const string& en : energies) {
        double end = sqrt(std::stod(en))*GeV;
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
    BinnedHistoPtr<string>  _npion;
    string ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(KLOE_2005_I655225);


}
