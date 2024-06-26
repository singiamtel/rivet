// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> K_1+/K_2+ phixs
  class BABAR_2008_I789278 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2008_I789278);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 10323);
      BP.addStable(-10323);
      BP.addStable(   325);
      BP.addStable(  -325);
      BP.addStable(   333);
      declare(BP, "BP");
      // histos
      for(unsigned int ix=0;ix<2;++ix) {
	book(_p[ix][0],1,1,1+ix);
        book(_p[ix][1],"TMP/wgt_"+toString(ix+1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 333,1}, { 10323,1}};
      static const map<PdgId,unsigned int> & mode1CC = { { 333,1}, {-10323,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 333,1}, {   325,1}};
      static const map<PdgId,unsigned int> & mode2CC = { { 333,1}, {  -325,1}};

      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
	int imode=-1;
       	if (BP.modeMatches(ix,2,mode1  ) ||
	    BP.modeMatches(ix,2,mode1CC) ) imode=0;
       	else if (BP.modeMatches(ix,2,mode2  ) ||
		 BP.modeMatches(ix,2,mode2CC) ) imode=1;
	else continue;
	// find phi decay products
	const Particle & phi = BP.decayProducts()[ix].at(333)[0];
	if(phi.children().size()!=2) continue;
	if(phi.children()[0].pid()!=-phi.children()[1].pid()) continue;
	if(phi.children()[0].abspid()!=321) continue;
	Particle Kp = phi.children()[0].pid()>0 ? phi.children()[0] : phi.children()[1];
      	// boost to B rest frame
	LorentzTransform boost =
       	  LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix].momentum().betaVec());
	FourMomentum pPhi = boost.transform(phi.momentum());
	FourMomentum pKp  = boost.transform( Kp.momentum());
	const LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pPhi.betaVec());
	pKp = boost2.transform(pKp);
	double cK = pKp.p3().unit().dot(pPhi.p3().unit());
	_p[imode][0]->fill(-(1.-5.*sqr(cK))/2.);
        _p[imode][1]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
        scale(_p[ix][0], 1./ *_p[ix][1]);
    }

    /// @}


    /// @name Histograms
    /// @{
   CounterPtr _p[2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2008_I789278);

}
