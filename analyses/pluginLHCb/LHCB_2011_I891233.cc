// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi production at 7 TeV
  class LHCB_2011_I891233 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2011_I891233);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<5;++ix) {
        book(_h_Jpsi[ix],{2.0,2.5,3.0,3.5,4.0,4.5});
        for(unsigned int iy=0;iy<5;++iy) {
          if(ix<4)
            book(_h_Jpsi[ix]->bin(iy+1),6+ix,1,iy+1);
          else
            book(_h_Jpsi[ix]->bin(iy+1),"TMP/JPsi_"+toString(iy),refData(6,1,iy+1));
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // J/psi
      for (const Particle& p : ufs.particles(Cuts::pid==443)) {
	// prompt/non-prompt
	bool nonPrompt = p.fromBottom();
        double absrap = p.absrap();
        double xp = p.perp();
	if(absrap<2. || absrap>4.5 ||  xp>14.) continue;
        _h_Jpsi[4]->fill(absrap,xp);
	if(nonPrompt) {
	  _h_Jpsi[1]->fill(absrap,xp);
	}
	else {
	  _h_Jpsi[0]->fill(absrap,xp);
	  _h_Jpsi[2]->fill(absrap,xp);
	  _h_Jpsi[3]->fill(absrap,xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      double factor = 0.5*crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<5;++ix) {
	scale(_h_Jpsi[ix],factor);
        divByGroupWidth(_h_Jpsi[ix]);
      }
      for(unsigned int ix=0;ix<_h_Jpsi[4]->numBins();++ix) {
	Estimate1DPtr tmp;
	book(tmp,10,1,ix+1);
	divide(_h_Jpsi[1]->bin(ix+1),_h_Jpsi[4]->bin(ix+1),tmp);
	tmp->scale(100.);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_Jpsi[5];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2011_I891233);

}
