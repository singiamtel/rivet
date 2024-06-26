// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief psi(2S) production at 7 and 13 TeV
  class LHCB_2019_I1748712 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2019_I1748712);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // histograms
      _is=-1;
      if (isCompatibleWithSqrtS(7000)) {
      	_is=1;
      }
      else if (isCompatibleWithSqrtS(13000)) {
      	_is=0;
      }
      else  {
      	throw Error("Invalid CMS energy for LHCB_2019_I1748712");
      }
      for(unsigned int ix=0;ix<3;++ix) {
        book(_h_psi[ix],{2.0,2.5,3.0,3.5,4.0,4.5});
        for(unsigned int iy=0;iy<5;++iy) {
          if(ix<2)
            book(_h_psi[ix]->bin(iy+1),1+ix+2*_is,1,iy+1);
          else
            book(_h_psi[ix]->bin(iy+1),"TMP/psi_"+toString(iy),refData(1+2*_is,1,iy+1));
        }
      }
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_pT[ix],5+2*_is,1,ix+1);
	book(_h_y [ix],6+2*_is,1,ix+1);
	if(_is==0) {
	  book(_h_pT_J[ix],"TMP/Jpsi_pT_"+toString(ix),refData(11,1,ix+1));
	  book(_h_pT_2[ix],"TMP/psi_pT_"+toString(ix),refData(11,1,ix+1));
	  book(_h_y_J [ix] ,"TMP/Jpsi_y_" +toString(ix),refData(12,1,ix+1));
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // J/psi
      for (const Particle& p : ufs.particles(Cuts::pid==443 || Cuts::pid==100443)) {
	// prompt/non-prompt
	bool nonPrompt = p.fromBottom();
        double absrap = p.absrap();
        double xp = p.perp();
	if(absrap<2. || absrap>4.5) continue;
	if     (_is==0 && (xp<2   || xp>20.)) continue;
	else if(_is==1 && (xp<3.5 || xp>14.)) continue;
	if (p.pid()==100443) {
	  _h_psi[nonPrompt]->fill(absrap,xp);
	  _h_psi[2        ]->fill(absrap,xp);
	  _h_pT[nonPrompt]->fill(xp);
	  _h_y [nonPrompt]->fill(absrap);
	  if(_h_pT_2[0]) _h_pT_2[nonPrompt]->fill(xp);
	}
	else if (_h_y_J[0]) {
	  _h_pT_J[nonPrompt]->fill(xp);
	  _h_y_J [nonPrompt]->fill(absrap);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      double factor = 0.5*crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<3;++ix) {
	scale(_h_psi[ix],factor);
        divByGroupWidth(_h_psi[ix]);
      }
      Estimate1DPtr tmp;
      for(unsigned int ix=0;ix<_h_psi[1]->numBins();++ix) {
	book(tmp,9+_is,1,1+ix);
	divide(_h_psi[1]->bin(ix+1),_h_psi[2]->bin(ix+1),tmp);
	tmp->scale(100.);
      }
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_pT[ix],factor);
	scale(_h_y [ix],factor);
      }
      // ratio psi(2s)/J/psi only at 13 TeV
      if(_is==0) {
	for(unsigned int ix=0;ix<2;++ix) {
	  scale(_h_pT_J[ix],factor);
	  scale(_h_pT_2[ix],factor);
	  scale(_h_y_J [ix],factor);
	  book(tmp,11,1,ix+1);
	  divide(_h_pT_2[ix],_h_pT_J[ix],tmp);
	  book(tmp,12,1,ix+1);
	  divide(_h_y[ix],_h_y_J[ix],tmp);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_psi[3];
    Histo1DPtr _h_pT[2],_h_pT_2[2],_h_y[2];
    Histo1DPtr _h_pT_J[2],_h_y_J[2];
    int _is;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2019_I1748712);

}
