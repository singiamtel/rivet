// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> Ds1(2460) Dbar
  class BABAR_2004_I656424 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2004_I656424);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511||
						Cuts::abspid==521);
      declare(ufs, "UFS");
      // histos
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles BB = apply<UnstableParticles>(event, "UFS").particles();
      for(const Particle & p : BB) {
	if(p.children().size()!=2) continue;
	int sign = p.pid()/p.abspid();
	int iD = -sign*(p.abspid()-100);
	Particle Ds1;
	if(p.children()[0].pid()==sign*20433 &&
	   p.children()[1].pid()==iD) {
	  Ds1=p.children()[0];
	}
	else if(p.children()[1].pid()==sign*20433 &&
		p.children()[0].pid()==iD) {
	  Ds1=p.children()[1];
	}
	else
	  continue;
	// find children of the Ds1 meson
	Particle Ds;
	if(Ds1.children().size()!=2) continue;
	if(Ds1.children()[0].abspid()==431 &&
	   Ds1.children()[1].abspid()==22)
	  Ds = Ds1.children()[0];
	else if (Ds1.children()[1].abspid()==431 &&
		 Ds1.children()[0].abspid()==22)
	  Ds = Ds1.children()[1];
	else
	  continue;
	// boost to rest frame
	LorentzTransform boostB  = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	FourMomentum pDs1  = boostB.transform(Ds1.momentum());
	FourMomentum pDs   = boostB.transform(Ds .momentum());
	LorentzTransform boostDs1 = LorentzTransform::mkFrameTransformFromBeta(pDs1.betaVec());
	pDs = boostDs1.transform(pDs);
	Vector3 axis = pDs1.p3().unit();
	double cTheta = axis.dot(pDs.p3().unit());
	_h->fill(cTheta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2004_I656424);

}
