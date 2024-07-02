// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> K* phi
  class BELLE_2005_I677873 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2005_I677873);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521);
      declare(ufs, "UFS");
      // histos
      for(unsigned int ix=0; ix<2; ++ix) {
        for(unsigned int iy=0; iy<3; ++iy) {
          book(_h[ix][iy], 1, 1+ix, 1+iy);
        }
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
      for (const Particle& B : ufs.particles()) {
        if (B.children().size()!=2) continue;
        Particle phi,Kstar;
        if (B.children()[0].pid()==PID::PHI &&
            isKstar(B.children()[1].abspid())) {
          phi = B.children()[0];
          Kstar = B.children()[1];
        }
        else if (B.children()[1].pid()==PID::PHI &&
                 isKstar(B.children()[0].abspid())) {
          phi = B.children()[1];
          Kstar = B.children()[0];
        }
        else {
          continue;
        }
        if (phi.children().size()!=2) continue;
        // find Kstar decay products
        Particle K;
        if (isK (Kstar.children()[0].abspid()) &&
            isPi(Kstar.children()[1].abspid())) {
          K = Kstar.children()[0];
        }
        else if (isK(Kstar.children()[1].abspid()) &&
                 isPi(Kstar.children()[0].abspid())) {
          K = Kstar.children()[1];
        }
        else {
          continue;
        }
        // find phi decay products
        Particle oDec;
        if (phi.children()[0].pid()==-321 && phi.children()[1].pid()== 321) {
          oDec = phi.children()[1];
        }
        else if (phi.children()[1].pid()==-321 && phi.children()[0].pid()== 321) {
          oDec = phi.children()[0];
        }
        else {
          continue;
        }
        // boost to B rest frame
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(B.mom().betaVec());
        FourMomentum pPhi = boost.transform(phi.mom());
        FourMomentum pKstar = boost.transform(Kstar.mom());
        FourMomentum pK     = boost.transform(K.mom());
        FourMomentum pOdec  = boost.transform(oDec.mom());
        // axes
        Vector3 axisX = pPhi.p3().unit();
        Vector3 axisY = (pK.p3()-axisX.dot(pK.p3())*axisX).unit();
        Vector3 axisZ = axisX.cross(axisY).unit();
        // kaon helicity angle
        LorentzTransform boostK = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
        double cosK = -axisX.dot(boostK.transform(pK).p3().unit());
        // transversality angles
        LorentzTransform boostL = LorentzTransform::mkFrameTransformFromBeta(pPhi.betaVec());
        Vector3 axisL = boostL.transform(pOdec).p3().unit();
        double cosL = axisL.dot(axisZ);
        double phiL = atan2(axisL.dot(axisY),axisL.dot(axisX));
        // fill hists
        int iloc = B.abspid()==511 ? 0 : 1;
        _h[iloc][0]->fill(cosK);
        _h[iloc][1]->fill(cosL);
        _h[iloc][2]->fill(phiL);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0; ix<2; ++ix) {
        normalize(_h[ix], 1.0, false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2005_I677873);

}
