// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief e+ e- > n nbar
  class SND_2014_I1321689 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2014_I1321689);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      for(unsigned int ix=0;ix<2;++ix)
        book(_nneutron[ix], 1+ix,1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if (fs.particles().size() != 2) vetoEvent;
      for (const Particle& p : fs.particles()) {
        if (p.abspid() != PID::NEUTRON) vetoEvent;
      }
      for(unsigned int ix=0;ix<2;++ix) _nneutron[ix]->fill(round(sqrtS()/MeV));
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for(unsigned int ix=0;ix<2;++ix)
        scale(_nneutron[ix], crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _nneutron[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_2014_I1321689);


}
