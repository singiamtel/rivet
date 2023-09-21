// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class DM1_1979_I132828 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DM1_1979_I132828);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // Book histograms
      book(_npion, 1, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (Ecm == "") {
        const size_t idx = _ecmAxis.index(sqrtS()/MeV);
        Ecm = idx? _npion->xEdges()[idx-1] : "OTHER";
      }
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=4) vetoEvent;
      for (const Particle& p : fs.particles()) {
        if(abs(p.pid())!=PID::PIPLUS)  vetoEvent;
      }
      _npion->fill(Ecm);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_npion, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _npion;
    YODA::Axis<double> _ecmAxis{963., 1008., 1024., 1055., 1088., 1108.};
    string Ecm = "";
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(DM1_1979_I132828);


}
