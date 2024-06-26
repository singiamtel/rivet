// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  B -> D(*) K- KS0
  class BELLE_2023_I2655951 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2023_I2655951);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 or Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 411);
      BB.addStable(-411);
      BB.addStable( 421);
      BB.addStable(-421);
      BB.addStable( 413);
      BB.addStable(-413);
      BB.addStable( 423);
      BB.addStable(-423);
      BB.addStable( 310);
      BB.addStable( 313);
      BB.addStable(-313);
      declare(BB, "BB");
      for(unsigned int ix=0;ix<4;++ix) {
        book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const double mDs = 1.96835;
      static const map<PdgId,unsigned int> modes  [4] = { { {-421,1}, { 321,1}, { 310,1} },
      							  { {-411,1}, { 321,1}, { 310,1} },
      							  { {-423,1}, { 321,1}, { 310,1} },
      							  { {-413,1}, { 321,1}, { 310,1} }};
      static const map<PdgId,unsigned int> modesCC[4] = { { { 421,1}, {-321,1}, { 310,1} },
      							  { { 411,1}, {-321,1}, { 310,1} },
      							  { { 423,1}, {-321,1}, { 310,1} },
      							  { { 413,1}, {-321,1}, { 310,1} }};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      for (unsigned int ix=0;ix<BB.decaying().size();++ix) {
      	int imode=-1;
      	for (unsigned int iy=0;iy<4;++iy) {
      	  if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,modes[iy])) {
      	    imode=iy;
      	    break;
      	  }
      	  else if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,modesCC[iy])) {
      	    imode=iy;
      	    break;
      	  }
      	}
      	if (imode<0) continue;
      	int sign = BB.decaying()[ix].pid()/BB.decaying()[ix].abspid();
        const Particle & Km = BB.decayProducts()[ix].at( sign*321)[0];
        const Particle & K0 = BB.decayProducts()[ix].at( 310     )[0];
        double mKK=(Km.momentum()+K0.momentum()).mass();
        if (abs(mKK-mDs)<0.02) continue;
        _h[imode]->fill(mKK);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2023_I2655951);

}
