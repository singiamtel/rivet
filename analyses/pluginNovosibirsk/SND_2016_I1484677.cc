// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+ e- > K+K-
  class SND_2016_I1484677 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2016_I1484677);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      for(unsigned int ix=0;ix<2;++ix) {
        book(_nkaon[ix], 1+ix, 1, 1);
        for (const string& en : _nkaon[ix].binning().edges<0>()) {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      if (_ecms[0].empty() && _ecms[1].empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
	if(abs(p.pid())!=PID::KPLUS) vetoEvent;
      }
      for(unsigned int ix=0;ix<2;++ix) _nkaon[ix]->fill(_ecms[ix]);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nkaon, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nkaon[2];
    string _ecms[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_2016_I1484677);


}
