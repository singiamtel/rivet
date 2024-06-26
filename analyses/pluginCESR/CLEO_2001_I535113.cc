// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> K* psi(2S)
  class CLEO_2001_I535113 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2001_I535113);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521);
      declare(ufs, "UFS");
      // histos
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h[ix],1,1,1+ix);
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
        Particle onium,Kstar;
        if (B.children()[0].pid()==100443 &&
           isKstar(B.children()[1].abspid())) {
          onium = B.children()[0];
          Kstar = B.children()[1];
        }
        else if (B.children()[1].pid()==100443 &&
          isKstar(B.children()[0].abspid())) {
          onium = B.children()[1];
          Kstar = B.children()[0];
        }
        else {
          continue;
        }
        // find Kstar decay products
        Particle K;
        if (isK (Kstar.children()[0].abspid()) &&
            isPi(Kstar.children()[1].abspid())) {
          K = Kstar.children()[0];
        }
        else if(isK (Kstar.children()[1].abspid()) &&
          isPi(Kstar.children()[0].abspid())) {
          K = Kstar.children()[1];
        }
        else {
          continue;
        }
        unsigned int imode=0;
        if (B.abspid()==511 && K.abspid()==321) imode=2;
        else if(B.abspid()==521) {
          if(K.abspid()==321) imode=1;
          else                imode=0;
        }
        else {
          continue;
        }
        // boost to B rest frame
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(B.mom().betaVec());
        FourMomentum pOnium = boost.transform(onium.mom());
        FourMomentum pKstar = boost.transform(Kstar.mom());
        FourMomentum pK     = boost.transform(K    .mom());
        // axes
        Vector3 axisX = pOnium.p3().unit();
        // kaon helicity angle
        LorentzTransform boostK = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
        double cosK = -axisX.dot(boostK.transform(pK).p3().unit());
        _h[imode]->fill(cosK);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2001_I535113);

}
