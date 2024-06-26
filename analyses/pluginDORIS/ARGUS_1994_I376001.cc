// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> J/psi K*
  class ARGUS_1994_I376001 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1994_I376001);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521);
      declare(ufs, "UFS");
      // histograms
      for (unsigned int iz=0; iz<2; ++iz) {
        book(_h[iz],1,1,1+iz);
      }
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
      for (const Particle & B : ufs.particles()) {
      	if (B.children().size()!=2) continue;
      	Particle jPsi,Kstar;
      	if (B.children()[0].pid()==443 && isKstar(B.children()[1].abspid())) {
      	  jPsi = B.children()[0];
      	  Kstar = B.children()[1];
      	}
      	else if (B.children()[1].pid()==443 && isKstar(B.children()[0].abspid())) {
      	  jPsi = B.children()[1];
      	  Kstar = B.children()[0];
      	}
      	else {
      	  continue;
        }
      	if (jPsi.children().size()!=2) continue;
      	// find Kstar decay products
      	Particle K;
      	if (isK (Kstar.children()[0].abspid()) && isPi(Kstar.children()[1].abspid())) {
      	  K = Kstar.children()[0];
      	}
      	else if (isK (Kstar.children()[1].abspid()) && isPi(Kstar.children()[0].abspid())) {
      	  K = Kstar.children()[1];
      	}
      	else {
      	  continue;
        }
      	// find jPsi decay products
      	Particle oDec;
      	if (jPsi.children()[0].pid()==-13 && jPsi.children()[1].pid()== 13) {
      	  oDec = jPsi.children()[1];
        }
        else if (jPsi.children()[1].pid()==-13 && jPsi.children()[0].pid()== 13) {
          oDec = jPsi.children()[0];
        }
        else if (jPsi.children()[0].pid()==-11 && jPsi.children()[1].pid()== 11) {
          oDec = jPsi.children()[1];
        }
        else if (jPsi.children()[1].pid()==-11 && jPsi.children()[0].pid()== 11) {
          oDec = jPsi.children()[0];
        }
        else {
          continue;
        }
      	// boost to B rest frame
      	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(B.mom().betaVec());
      	FourMomentum pJPsi = boost.transform(jPsi.mom());
      	FourMomentum pKstar = boost.transform(Kstar.mom());
      	FourMomentum pK     = boost.transform(K    .mom());
      	FourMomentum pOdec  = boost.transform(oDec .mom());
      	// axes
      	Vector3 axisX = pJPsi.p3().unit();
      	Vector3 axisY = (pK.p3()-axisX.dot(pK.p3())*axisX).unit();
      	Vector3 axisZ = axisX.cross(axisY).unit();
      	// kaon helicity angle
      	LorentzTransform boostK = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
      	double cosK = -axisX.dot(boostK.transform(pK).p3().unit());
      	// transversality angles
      	LorentzTransform boostL = LorentzTransform::mkFrameTransformFromBeta(pJPsi.betaVec());
      	Vector3 axisL = boostL.transform(pOdec).p3().unit();
      	double cosL = axisL.dot(axisZ);
      	// fill hists
       	_h[0]->fill(cosL);
       	_h[1]->fill(cosK);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1994_I376001);

}
