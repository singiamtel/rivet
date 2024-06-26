// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- > pi0 pi0 gamma
  class KLOE_2007_I725483 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(KLOE_2007_I725483);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // get the CMS point
      unsigned int ix=0;
      if (isCompatibleWithSqrtS(1.01955))     ix=1;
      else if(isCompatibleWithSqrtS(1.01965)) ix=2;
      else if(isCompatibleWithSqrtS(1.01985)) ix=3;
      else if(isCompatibleWithSqrtS(1.01995)) ix=4;
      else throw Error("Unexpected sqrtS !");
      for (unsigned int iy=0; iy<2; ++iy) {
        book(_h[iy],1,ix,iy+1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      Particles pi0,gamma;
      map<long,Particles> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        if (p.pid()==PID::PI0)          pi0.push_back(p);
        else if (p.pid()==PID::GAMMA) gamma.push_back(p);
        ++ntotal;
      }
      // three particles (pi0 pi0 gamma)
      if (ntotal!=3 || pi0.size()!=2 || gamma.size()!=1) vetoEvent;
      _h[0]->fill((pi0[0].mom()+  pi0[1].mom()).mass()/MeV);
      _h[1]->fill((pi0[0].mom()+gamma[0].mom()).mass()/MeV);
      _h[1]->fill((pi0[1].mom()+gamma[0].mom()).mass()/MeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(KLOE_2007_I725483);

}
