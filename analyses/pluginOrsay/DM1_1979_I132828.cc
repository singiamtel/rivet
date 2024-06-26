// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- -> pi+pi+pi-pi-
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
      // extract the centre-of-mass energies
      vector<string> energies = _npion.binning().edges<0>();
      if     (sqrtS()>=0.963 && sqrtS()<=1.002)  _ecms=energies[0];
      else if(sqrtS()>=1.008 && sqrtS()<=1.024)  _ecms=energies[1];
      else if(isCompatibleWithSqrtS(std::stod(energies[2])*MeV)) _ecms=energies[2];
      else if(isCompatibleWithSqrtS(std::stod(energies[3])*MeV)) _ecms=energies[3];
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=4) vetoEvent;
      for (const Particle& p : fs.particles()) {
        if(abs(p.pid())!=PID::PIPLUS)  vetoEvent;
      }
      _npion->fill(_ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_npion, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _npion;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DM1_1979_I132828);


}
