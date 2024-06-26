// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- to Upsilon to mu+mu-
  class CRYSTAL_BALL_1992_I306832 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CRYSTAL_BALL_1992_I306832);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
        book(_c[ix], 1+ix, 1, 1);
      _eCent = getOption<string>("ECENT", std::to_string(sqrtS()/MeV));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over FS particles
      Particles fs = apply<FinalState>(event,"FS").particles();
      Particles mm,mp;
      for (const Particle& p : fs) {
        if (p.abspid()==PID::MUON) {
          if (p.pid()>0) mm.push_back(p);
          else           mp.push_back(p);
        }
        else if(p.pid()!=PID::GAMMA) {
          vetoEvent;
        }
      }
      if (mm.size()==1 && mp.size()==1) {
        _c[0]->fill(_eCent);
        _c[1]->fill(_eCent);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      for (unsigned int ix=0; ix<2; ++ix)
        scale(_c[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c[2];
    string _eCent;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CRYSTAL_BALL_1992_I306832);

}
