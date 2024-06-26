// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> p pbar K+
  class BABAR_2005_I686573 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2005_I686573);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 10441);
      BP.addStable( 10443);
      BP.addStable( 20443);
      BP.addStable(100443);
      BP.addStable(   441);
      BP.addStable(   445);
      declare(BP, "BP");
      // histos
      book(_h_pp,1,1,1);
      for(unsigned int ix=0;ix<2;++ix)
	book(_h_pK[ix],2,1,1+ix);
      book(_c,"TMP/nB");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 2212,1}, {-2212,1}, { 321,1} };
      static const map<PdgId,unsigned int> & mode1CC = { { 2212,1}, {-2212,1}, {-321,1} };
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
       	int sign=1;
	if(BP.modeMatches(ix,3,mode1))        sign = 1;
	else if(BP.modeMatches(ix,3,mode1CC)) sign =-1;
	else continue;
	_c->fill();
	const Particle & pp   = BP.decayProducts()[ix].at( sign*2212)[0];
	const Particle & pbar = BP.decayProducts()[ix].at(-sign*2212)[0];
	const Particle & Kp   = BP.decayProducts()[ix].at( sign*321 )[0];
	_h_pp->fill((pp.momentum()+pbar.momentum()).mass());
	double mpK    = (Kp.momentum()+pp  .momentum()).mass();
	double mpbarK = (Kp.momentum()+pbar.momentum()).mass();
	if(mpK>mpbarK) _h_pK[0]->fill(mpK   );
	else           _h_pK[1]->fill(mpbarK);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pp,1./ *_c);
      for(unsigned int ix=0;ix<2;++ix)
       	scale(_h_pK[ix],1./ *_c);
      Histo1DPtr tmp;
      book(tmp,2,1,3);
      *tmp = *_h_pK[0] - *_h_pK[1];
      tmp->setPath("/"+name()+"/"+mkAxisCode(2,1,3));
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pp,_h_pK[2];
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2005_I686573);

}
