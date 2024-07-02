// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+ e- > n nbar
  class SND_2022_I2102082 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2022_I2102082);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_nneutron, 1, 1, 1);

      for (const string& en : _nneutron.binning().edges<0>()) {
        double end = std::stod(en)*MeV;
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
      if (fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
        if(abs(p.pid())!=PID::NEUTRON) vetoEvent;
      }
      _nneutron->fill(_ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nneutron, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nneutron;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_2022_I2102082);

}
