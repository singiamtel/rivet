// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> rho+ rho-
  class BABAR_2004_I632736 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2004_I632736);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable( 213);
      B0.addStable(-213);
      declare(B0, "B0");
      // histos
      book(_p[0],1,1,1);
      book(_p[1],"TMP/wgt");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 213,1}, {-213,1} };
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
       	int sign = B0.decaying()[ix].pid()/B0.decaying()[ix].abspid();
      	if(!B0.modeMatches(ix,2,mode)) continue;
	const Particle & rhop = B0.decayProducts()[ix].at( sign*213)[0];
	//const Particle & rhom = B0.decayProducts()[ix].at(-sign*213)[0];
	// compute the helicity angle
	if(rhop.children().size()!=2) continue;
	Particle pip;
	if(rhop.children()[0].pid()==sign*211 &&
	   rhop.children()[1].pid()==111)
	  pip = rhop.children()[0];
	else if(rhop.children()[1].pid()==sign*211 &&
		rhop.children()[0].pid()==111)
	  pip = rhop.children()[1];
	  else
	    continue;
	LorentzTransform boost1 =
	  LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix].momentum().betaVec());
	FourMomentum prho = boost1.transform(rhop.momentum());
	FourMomentum ppi  = boost1.transform(pip);
	LorentzTransform boost2 =
	  LorentzTransform::mkFrameTransformFromBeta(prho.betaVec());
	ppi = boost2.transform(ppi);
	double cTheta = ppi.p3().unit().dot(prho.p3().unit());
	_p[0]->fill(-0.5*(1.-5.*sqr(cTheta)));
        _p[1]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_p[0],1./ *_p[1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _p[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2004_I632736);

}
