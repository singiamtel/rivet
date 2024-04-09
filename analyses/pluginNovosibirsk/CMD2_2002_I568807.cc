// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class CMD2_2002_I568807 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD2_2002_I568807);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_npion, 1,1,2);
      vector<string> energies({"610.5", "620.5", "630.5", "640.51", "650.49", "660.5", "670.5", "680.59", "690.43", "700.52",
          "710.47", "720.25", "730.24", "740.2", "750.28", "760.18", "764.17", "770.11", "774.38", "778.17", "780.17", "782.23",
          "784.24", "786.04", "790.1", "794.14", "800.02", "810.14", "820.02", "829.97", "839.1", "849.24", "859.6", "869.5",
          "879.84", "889.72", "900.04", "910.02", "919.56", "930.11", "942.19", "951.84", "961.52"});
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


  RIVET_DECLARE_PLUGIN(CMD2_2002_I568807);


}
