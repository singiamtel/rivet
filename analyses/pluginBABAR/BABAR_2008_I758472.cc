// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Longitudinal Polarization in B -> K*0 K*0
  class BABAR_2008_I758472 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2008_I758472);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 313);
      BB.addStable(-313);
      declare(BB, "BB");
      // histos
      book(_p[0],1,1,1);
      book(_p[1],"TMP/wgt");
    }

    bool isK(int id) {
      return id==321 || id==311 || id==310 || id==130; 
    }

    bool isPi(int id) {
      return id==211 || id==111;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { {-313,1}, { 313,1}};

      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
	if (!BB.modeMatches(ix,2,mode1)) continue;
	int sign = BB.decaying()[ix].pid()/BB.decaying()[ix].abspid();
       	const Particle & Kstar = BB.decayProducts()[ix].at( sign*313)[0];
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
	  LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix]. momentum().betaVec());
        FourMomentum pKstar  = boost.transform(Kstar.momentum());
        FourMomentum pK = boost.transform(KK   .momentum());
        const LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
        pK = boost2.transform(pK);
        FourMomentum pB = boost2.transform(boost.transform(BB.decaying()[ix].momentum()));
        double cosK = -pB.p3().unit().dot(pK.p3().unit());
        _p[0]->fill(0.5*(5.*sqr(cosK)-1));
        _p[1]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_p[0], 1./ *_p[1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _p[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2008_I758472);

}
