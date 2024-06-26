// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> Ds* D*
  class BABAR_2003_I613283 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2003_I613283);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable( 413);
      B0.addStable(-413);
      B0.addStable( 433);
      B0.addStable(-433);
      declare(B0, "B0");
      // histos
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],2,1+ix,1);
      for(unsigned int ix=0;ix<2;++ix)
        book(_p[ix],"TMP/p_"+toString(ix+1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 433,1}, {-413,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-433,1}, { 413,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
	int sign = 1;
       	if      ( B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,2,mode  ) ) sign= 1.;
       	else if ( B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,2,modeCC) ) sign=-1.;
	else continue;
	// find D* decay products
	const Particle & Dstar = B0.decayProducts()[ix].at(-sign*413)[0];
	if(Dstar.children().size()!=2) continue;
	Particle pim;
	if(Dstar.children()[0].pid()==-sign*211 &&
	   Dstar.children()[1].pid()==-sign*421) {
	  pim = Dstar.children()[0];
	}
	else if(Dstar.children()[1].pid()==-sign*211 &&
		Dstar.children()[0].pid()==-sign*421) {
	  pim = Dstar.children()[1];
	}
	else
	  continue;
	// find Ds* decay products
	const Particle & DSstar = B0.decayProducts()[ix].at(sign*433)[0];
	Particle gamma;
	if(DSstar.children()[0].pid()   == 22 &&
	   DSstar.children()[1].pid()==sign*431) {
	  gamma= DSstar.children()[0];
	}
	else if(DSstar.children()[1].pid()   == 22 &&
	 	DSstar.children()[0].pid()==sign*431) {
	  gamma = DSstar.children()[1];
	}
	else
	  continue;
      	// boost to B rest frame
	LorentzTransform boost =
       	  LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix].momentum().betaVec());
	FourMomentum pDSstar = boost.transform(DSstar.momentum());
	FourMomentum pDstar  = boost.transform( Dstar.momentum());
	FourMomentum pGamma = boost.transform(gamma.momentum());
	FourMomentum pPim   = boost.transform(pim.momentum());
	const LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDstar.betaVec());
	pPim = boost2.transform(pPim);
	double cPi = pPim.p3().unit().dot(pDSstar.p3().unit());
	const LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pDSstar.betaVec());
	 pGamma = boost3.transform(pGamma);
	double cGamma = pGamma.p3().unit().dot(pDstar.p3().unit());
	_p[0]->fill(-100.*(1.-5.*sqr(cPi))/2.);
        _p[1]->fill();
	_h[0]->fill(cPi);
	_h[1]->fill(cGamma);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix],1.,false);
      Estimate0DPtr tmp;
      book(tmp,1,1,1);
      divide(*_p[0],*_p[1],tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    CounterPtr _p[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2003_I613283);

}
