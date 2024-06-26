// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> K* J/psi, psi(2S) and chi_c1
  class BABAR_2007_I748372 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2007_I748372);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 or
						Cuts::abspid==521);
      declare(ufs, "UFS");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<3;++iy)
	  for(unsigned int iz=0;iz<3;++iz)
	    book(_h[ix][iy][iz],1+ix,1+iy,1+iz);
    }

    bool isOnium(int pid) const {
      return pid==443 || pid==100443 || pid==20443;
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
	Particle onium,Kstar;
	if(isOnium(B.children()[0].abspid()) &&
	   isKstar(B.children()[1].abspid())) {
	  onium = B.children()[0];
	  Kstar = B.children()[1];
	}
	else if(isOnium(B.children()[1].abspid()) &&
		isKstar(B.children()[0].abspid())) {
	  onium = B.children()[1];
	  Kstar = B.children()[0];
	}
	else
	  continue;
	if(onium.children().size()!=2) continue;
	// find Kstar decay products
	Particle K;
	if(isK (Kstar.children()[0].abspid()) &&
	   isPi(Kstar.children()[1].abspid())) {
	  K = Kstar.children()[0];
	}
	else if(isK (Kstar.children()[1].abspid()) &&
		isPi(Kstar.children()[0].abspid())) {
	  K = Kstar.children()[1];
	}
	else
	  continue;
	int iK=-1;
	if(B.abspid()==511) {
	  if(K.abspid()==321) iK=0;
	  else continue;
	}
	else {
	  if(K.abspid()==321) iK=2;
	  else                iK=1;
	}
	// find onium decay products
	Particle oDec;
	if(onium.pid()==20443) {
	  if(onium.children()[0].pid()==443 &&
	     onium.children()[1].pid()==22) {
	    oDec = onium.children()[0];
	  }
	  else if(onium.children()[1].pid()==443 &&
		  onium.children()[0].pid()==22) {
	    oDec = onium.children()[1];
	  }
	  else
	    continue;
	}
	else {
	  if(onium.children()[0].pid()==-13 &&
	     onium.children()[1].pid()== 13) {
	    oDec = onium.children()[1];
	  }
	  else if(onium.children()[1].pid()==-13 &&
		  onium.children()[0].pid()== 13) {
	    oDec = onium.children()[0];
	  }
	  else if(onium.children()[0].pid()==-11 &&
	     onium.children()[1].pid()== 11) {
	    oDec = onium.children()[1];
	  }
	  else if(onium.children()[1].pid()==-11 &&
		  onium.children()[0].pid()== 11) {
	    oDec = onium.children()[0];
	  }
	  else
	    continue;
	}
	// boost to B rest frame
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(B.momentum().betaVec());
	FourMomentum pOnium = boost.transform(onium.momentum());
	FourMomentum pKstar = boost.transform(Kstar.momentum());
	FourMomentum pK     = boost.transform(K    .momentum());
	FourMomentum pOdec  = boost.transform(oDec .momentum());
	// axes
	Vector3 axisX = pOnium.p3().unit();
	Vector3 axisY = (pK.p3()-axisX.dot(pK.p3())*axisX).unit();
	Vector3 axisZ = axisX.cross(axisY).unit();
	// kaon helicity angle
	LorentzTransform boostK = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
	double cosK = -axisX.dot(boostK.transform(pK).p3().unit());
	// transversality angles
	LorentzTransform boostL = LorentzTransform::mkFrameTransformFromBeta(pOnium.betaVec());
	Vector3 axisL = boostL.transform(pOdec).p3().unit();
	double cosL = axisL.dot(axisZ);
	double phiL = atan2(axisL.dot(axisY),axisL.dot(axisX));
	if(phiL<0.) phiL +=2.*M_PI;
	// fill hists
	int iloc = onium.pid()==443 ? 0 : (onium.pid()==100443 ? 1 : 2);
	_h[0][iloc][iK]->fill(cosK);
	_h[1][iloc][iK]->fill(cosL);
	_h[2][iloc][iK]->fill(phiL);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<3;++iy)
	  for(unsigned int iz=0;iz<3;++iz)
	    normalize(_h[ix][iy][iz],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3][3][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2007_I748372);

}
