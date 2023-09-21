// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> K* l+ l-
  class BELLE_2019_I1728173 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2019_I1728173);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 or
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(   443);
      BB.addStable(100443);
      BB.addStable( 313);
      BB.addStable( 323);
      BB.addStable(-313);
      BB.addStable(-323);
      declare(BB, "BB");
      // book histograms
      for(unsigned int ix=0;ix<4;++ix)
      	for(unsigned int iy=0;iy<3;++iy) {
      	  book(_h_br[ix][iy],"TMP/h_br_"+toString(ix+1)+"_"+toString(iy+1),refData(1,1+iy,1));
      	}
      for(unsigned int ix=0;ix<2;++ix)
      	book(_c[ix],"TMP/nB_"+toString(ix+1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 323,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-323,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 313,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {-313,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode3   = { { 323,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode3CC = { {-323,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode4   = { { 313,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode4CC = { {-313,1},{ 11,1}, {-11,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
      	if(BB.decaying()[ix].abspid()==521) _c[0]->fill();
      	else                                _c[1]->fill();
      	int imode=0;
      	if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode1)) ||
	    (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode1CC)))       imode=0;
	else if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode2)) ||
		 (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode2CC)))  imode=1;
      	else if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode3)) ||
		 (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode3CC)))  imode=2;
      	else if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode4)) ||
		 (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode4CC)))  imode=3;
      	else continue;
      	int il = imode<2 ? 13 : 11;
      	const Particle & lp = BB.decayProducts()[ix].at(-il)[0];
      	const Particle & lm = BB.decayProducts()[ix].at( il)[0];
      	double qq = (lp.momentum()+lm.momentum()).mass2();
      	for(unsigned int iy=0;iy<3;++iy) {
      	  _h_br[imode][iy]->fill(qq);
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // ratio of lifetimes
      double rLife = 1.078;
      // normalize BR plots
      for(unsigned int ix=0;ix<4;++ix) {
      	for(unsigned int iy=0;iy<3;++iy) {
      	  if(ix%2==0) {
      	     scale(_h_br[ix][iy],1e7/ *_c[0]);
      	  }
      	  else {
      	     scale(_h_br[ix][iy],1e7*rLife/ *_c[1]);
      	  }
      	}
      }
      // RK*
      for (unsigned int ix=0;ix<3;++ix) {
      	Estimate1DPtr RK;
      	book(RK,1,1+ix,3);
      	divide(_h_br[0][ix],_h_br[2][ix],RK);
      	book(RK,1,1+ix,2);
      	divide(_h_br[1][ix],_h_br[3][ix],RK);
      	book(RK,1,1+ix,1);
      	for (size_t ibin=1; ibin < _h_br[1][ix]->numBins()+1; ++ibin) {
      	  double num     = _h_br[0][ix]->bin(ibin).sumW()   +_h_br[1][ix]->bin(ibin).sumW();
      	  double numErr2 = sqr(_h_br[0][ix]->bin(ibin).errW())+sqr(_h_br[1][ix]->bin(ibin).errW());
      	  double den     = _h_br[2][ix]->bin(ibin).sumW()   +_h_br[3][ix]->bin(ibin).sumW();
      	  double denErr2 = sqr(_h_br[2][ix]->bin(ibin).errW())+sqr(_h_br[3][ix]->bin(ibin).errW());
      	  double val(0.),err(0.);
      	  if(num>0. && den>0.) {
      	    val = num/den;
      	    err = val*(numErr2/sqr(num)+denErr2/sqr(den));
      	  }
      	  RK->bin(ibin).set(val, err);
      	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c[2];
    Histo1DPtr _h_br[4][3];
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BELLE_2019_I1728173);

}
