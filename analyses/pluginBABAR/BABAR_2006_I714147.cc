// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B- -> D0 K*-
  class BABAR_2006_I714147 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I714147);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      // hists
      book(_h,1,1,1);
    }

    bool isKstar(int pid) const {
      return pid==313 || pid==323;
    }

    bool isK(int pid) const {
      return pid==130 || pid==310 || pid==311 || pid==321;
    }

    bool isPi(int pid) const {
      return pid==211 || pid==111;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      UnstableParticles ufs = apply<UnstableParticles>(event, "UFS");
      for(const Particle & B : ufs.particles()) {
	if(B.children().size()!=2) continue;
	int sign = B.pid()/B.abspid();
	Particle D,Kstar;
	if(B.children()[0].pid()==-sign*421 &&
	   isKstar(B.children()[1].abspid())) {
	  D     = B.children()[0];
	  Kstar = B.children()[1];
	}
	else if(B.children()[1].pid()==-sign*421 &&
		isKstar(B.children()[0].abspid())) {
	  D     = B.children()[1];
	  Kstar = B.children()[0];
	}
	else
	  continue;
	// find Kstar decay products
	Particle pi;
	if(isK (Kstar.children()[0].abspid()) &&
	   isPi(Kstar.children()[1].abspid())) {
	  pi = Kstar.children()[1];
	}
	else if(isK (Kstar.children()[1].abspid()) &&
		isPi(Kstar.children()[0].abspid())) {
	  pi = Kstar.children()[0];
	}
	else
	  continue;
	// boost to B rest frame
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(B.momentum().betaVec());
	FourMomentum pD     = boost.transform(D    .momentum());
	FourMomentum pKstar = boost.transform(Kstar.momentum());
	FourMomentum pPi     = boost.transform(pi   .momentum());
	// pion helicity angle
	Vector3 axisX = pD.p3().unit();
	LorentzTransform boostK = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
	double cosPi = -axisX.dot(boostK.transform(pPi).p3().unit());
	_h->fill(cosPi);
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


  RIVET_DECLARE_PLUGIN(BABAR_2006_I714147);

}
