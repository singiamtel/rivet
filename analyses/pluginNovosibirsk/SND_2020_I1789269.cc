// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- -> pi+pi-
  class SND_2020_I1789269 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2020_I1789269);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_npion, 1,1,1);
      vector<string> energies({"525.1", "544.0", "565.2", "585.0", "604.8", "624.8", "644.6", "664.5", "684.4", "704.2", "724.1",
          "739.1", "743.8", "747.7", "751.7", "755.7", "759.6", "763.6", "767.8", "771.6", "775.7", "778.6", "780.7", "782.0",
          "782.9", "783.7", "784.7", "786.7", "789.5", "793.9", "797.7", "804.0", "821.8", "843.4", "862.7", "883.2"});
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


  RIVET_DECLARE_PLUGIN(SND_2020_I1789269);

}
