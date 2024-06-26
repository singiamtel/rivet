// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  B -> p pbar meson
  class BABAR_2007_I755548 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2007_I755548);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 ||
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 310);
      BB.addStable( 313);
      BB.addStable( 323);
      BB.addStable(-313);
      BB.addStable(-323);
      BB.addStable(   441);
      BB.addStable(   443);
      BB.addStable(100443);
      BB.addStable( 10441);
      BB.addStable( 20443);
      BB.addStable(   445);
      declare(BB, "BB");
      // histograms
      for(unsigned int ix=0;ix<4;++ix) {
	book(_h_mass[ix],1,1,1+ix);
	book(_c[ix],"TMP/c_"+toString(ix+1));
	for(unsigned int iy=0;iy<2;++iy)
	  if(!(ix==0&&iy==1)) book(_h_mass2[iy][ix],2,1+iy,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 2212,1}, {-2212,1}, { 310,1} };
      static const map<PdgId,unsigned int> & mode1CC = { { 2212,1}, {-2212,1}, { 310,1} };
      static const map<PdgId,unsigned int> & mode2   = { { 2212,1}, {-2212,1}, { 321,1} };
      static const map<PdgId,unsigned int> & mode2CC = { { 2212,1}, {-2212,1}, {-321,1} };
      static const map<PdgId,unsigned int> & mode3   = { { 2212,1}, {-2212,1}, { 313,1} };
      static const map<PdgId,unsigned int> & mode3CC = { { 2212,1}, {-2212,1}, {-313,1} };
      static const map<PdgId,unsigned int> & mode4   = { { 2212,1}, {-2212,1}, { 323,1} };
      static const map<PdgId,unsigned int> & mode4CC = { { 2212,1}, {-2212,1}, {-323,1} };
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
	int sign = BB.decaying()[ix].pid()/BB.decaying()[ix].abspid();
	int imode=-1, iK;
	if((sign>0 &&BB.modeMatches(ix,3,mode1)) ||
	   (sign<0 &&BB.modeMatches(ix,3,mode1CC))) {
	  imode = 0;
	  iK = 310;
	}
	else if((sign>0 &&BB.modeMatches(ix,3,mode2)) ||
		(sign<0 &&BB.modeMatches(ix,3,mode2CC))) {
	  imode = 1;
	  iK = sign*321;
	}
	else if((sign>0 &&BB.modeMatches(ix,3,mode3)) ||
		(sign<0 &&BB.modeMatches(ix,3,mode3CC))) {
	  imode = 2;
	  iK = sign*313;
	}
	else if((sign>0 &&BB.modeMatches(ix,3,mode4)) ||
		(sign<0 &&BB.modeMatches(ix,3,mode4CC))) {
	  imode = 3;
	  iK = sign*323;
	}
	else continue;
	_c[imode]->fill();
	const Particle & pp    = BB.decayProducts()[ix].at( sign*2212)[0];
	const Particle & pbar  = BB.decayProducts()[ix].at(-sign*2212)[0];
	const Particle & meson = BB.decayProducts()[ix].at(        iK)[0];
	FourMomentum pbaryon = pp.momentum()+pbar.momentum();
	double mass = pbaryon.mass();
	_h_mass[imode]->fill(mass);
	double mph[2] = {(pp  .momentum()+meson.momentum()).mass(),
			 (pbar.momentum()+meson.momentum()).mass()};
	if(imode==0) {
	  _h_mass2[0][imode]->fill(max(mph[0],mph[1]));
	}
	else {
	  if(mph[0]>mph[1])
	    _h_mass2[0][imode]->fill(mph[0]);
	  else 
	    _h_mass2[1][imode]->fill(mph[1]);
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix) {
	normalize(_h_mass[ix],1.,false);
	for(unsigned int iy=0;iy<2;++iy)
	  if(_h_mass2[iy][ix]) scale(_h_mass2[iy][ix],1./ *_c[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[4],_h_mass2[2][4];
    CounterPtr _c[4];
    /// @}
  };


  RIVET_DECLARE_PLUGIN(BABAR_2007_I755548);

}
