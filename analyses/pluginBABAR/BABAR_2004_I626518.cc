// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B- -> D*0 K*-
  class BABAR_2004_I626518 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2004_I626518);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 423);
      BP.addStable(-423);
      BP.addStable( 323);
      BP.addStable(-323);
      declare(BP, "BP");
      // histos
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[ix][iy],2,1+ix,1+iy);
      for(unsigned int ix=0;ix<2;++ix)
        book(_p[ix],"TMP/p_"+toString(ix+1));
    }

    bool isK(int id) {
      return id==321 || id==311 || id==310 || id==130; 
    }

    bool isPi(int id) {
      return id==211 || id==111;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {-423,1}, { 323,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 423,1}, {-323,1}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
	int sign = 1;
       	if      ( BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,2,mode  ) ) sign= 1.;
       	else if ( BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,2,modeCC) ) sign=-1.;
	else continue;
	// find D* decay products
	const Particle & Dstar = BP.decayProducts()[ix].at(-sign*423)[0];
	if(Dstar.children().size()!=2) continue;
	Particle D0;
	unsigned int imode=0;
	if(Dstar.children()[0].pid()== 111 &&
	   Dstar.children()[1].pid()==-sign*421) {
	  D0 = Dstar.children()[1];
	  imode=0;
	}
	else if(Dstar.children()[1].pid()== 111 &&
		Dstar.children()[0].pid()==-sign*421) {
	  imode=0;
	  D0 = Dstar.children()[0];
	}
	else if(Dstar.children()[0].pid()== 22 &&
		Dstar.children()[1].pid()==-sign*421) {
	  D0 = Dstar.children()[1];
	  imode=1;
	}
	else if(Dstar.children()[1].pid()== 22 &&
		Dstar.children()[0].pid()==-sign*421) {
	  imode=1;
	  D0 = Dstar.children()[0];
	}
	else
	  continue;
	// K* decay products
	const Particle & Kstar = BP.decayProducts()[ix].at( sign*323)[0];
	if(Kstar.children().size()!=2) continue;
       	Particle KK;
      	if(isK (Kstar.children()[0].abspid()) &&
      	   isPi(Kstar.children()[1].abspid()))
      	  KK = Kstar.children()[0];
      	else if(isK (Kstar.children()[1].abspid()) &&
      		isPi(Kstar.children()[0].abspid()))
      	  KK = Kstar.children()[1];
      	else continue;
      	// boost to B rest frame
	LorentzTransform boost =
       	  LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix]. momentum().betaVec());
	FourMomentum pKstar  = boost.transform(Kstar.momentum());
	FourMomentum pDstar  = boost.transform(Dstar.momentum());
	FourMomentum pK = boost.transform(KK   .momentum());
	FourMomentum pD = boost.transform(D0   .momentum());
	const LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDstar.betaVec());
	pD = boost2.transform(pD);
	double cD = pD.p3().unit().dot(pDstar.p3().unit());
	const LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
	pK = boost3.transform(pK);
	double cK = pK.p3().unit().dot(pKstar.p3().unit());
	if(imode==0) {
          _p[0]->fill(-0.5*(1-5*sqr(cK)));
          _p[1]->fill();
        }
	_h[0][imode]->fill(cD);
	_h[1][imode]->fill(cK);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	normalize(_h[ix][iy],1.,false);
      Estimate0DPtr tmp;
      book(tmp,1,1,1);
      divide(*_p[0],*_p[1],tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][2];
    CounterPtr _p[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2004_I626518);

}
