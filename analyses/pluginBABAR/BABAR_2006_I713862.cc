// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Bbar0 -> D*+ omega pi-
  class BABAR_2006_I713862 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I713862);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(PID::PI0);
      BB.addStable( 413);
      BB.addStable(-413);
      BB.addStable(223);
      declare(BB, "BB");
      // histos
      for(unsigned int ix=0;ix<2;++ix)
      	book(_h[ix],1+ix,1,1);
      book(_p,3,1,1);
      book(_c,"TMP/nB");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { -413,1}, { 211,1}, { 223,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {  413,1}, {-211,1}, { 223,1}};
      // loop over particles
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
	_c->fill();
	int sign = BB.decaying()[ix].pid()/BB.decaying()[ix].abspid();
	if ( (sign== 1 && BB.modeMatches(ix,3,mode1) ) ||
	     (sign==-1 && BB.modeMatches(ix,3,mode1CC) ) ) {
	  const Particle & omega = BB.decayProducts()[ix].at(      223)[0];
	  const Particle & pim   = BB.decayProducts()[ix].at( sign*211)[0];
	  const Particle & Dstar = BB.decayProducts()[ix].at(-sign*413)[0];
	  double mX = (omega.momentum()+pim.momentum()).mass();
	  _h[0]->fill(sqr(mX));
	  if(Dstar.children().size()!=2) continue;
	  // boost to B rest frame
	  LorentzTransform boost1 =
	    LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix].momentum().betaVec());
	  FourMomentum pDpi = boost1.transform(Dstar.momentum()+pim.momentum());
	  LorentzTransform boost2 =
	    LorentzTransform::mkFrameTransformFromBeta(pDpi.betaVec());
	  FourMomentum pDstar = boost2.transform(boost1.transform(Dstar.momentum()));
	  double c1 = pDstar.p3().unit().dot(pDpi.p3().unit());
	  if(c1<0.5) _h[1]->fill(pDpi.mass());
	  Particle D0;
	  if(Dstar.children()[0].pid()   ==-sign*211 &&
	     Dstar.children()[1].abspid()==-sign*421)
	    D0 = Dstar.children()[1];
	  else if(Dstar.children()[1].pid()     ==-sign*211 &&
	  	    Dstar.children()[0].abspid()==-sign*421)
	    D0 = Dstar.children()[0];
	  else
	    continue;
	  FourMomentum pD0    = boost2.transform(boost1.transform(D0   .momentum()));
	  LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pDstar.betaVec());
	  pD0 = boost3.transform(pD0);
	  double cTheta = pD0.p3().unit().dot(pDstar.p3().unit());
	  _p->fill(mX,.5*(5.*sqr(cTheta)-1.));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // first hist is differential BR/ BR for leptons
      double br = .1033;
      scale(_h[0], 1./br/ *_c);
      normalize(_h[1],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    Profile1DPtr _p;
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2006_I713862);

}
