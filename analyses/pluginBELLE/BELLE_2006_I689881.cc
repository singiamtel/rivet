// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> chi_c1 K*
  class BELLE_2006_I689881 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2006_I689881);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      // histograms
      book(_h,1,1,1);
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
        if (B.children().size()!=2) continue;
        Particle chic,Kstar;
        if (B.children()[0].pid()==20443 && B.children()[1].abspid()==313) {
          chic = B.children()[0];
          Kstar = B.children()[1];
        }
        else if (B.children()[1].pid()==20443 && B.children()[0].abspid()==313) {
          chic = B.children()[1];
          Kstar = B.children()[0];
        }
        else {
          continue;
        }
        // find Kstar decay products
        Particle K;
        if (isK (Kstar.children()[0].abspid()) && isPi(Kstar.children()[1].abspid())) {
          K = Kstar.children()[0];
        }
        else if (isK (Kstar.children()[1].abspid()) &&isPi(Kstar.children()[0].abspid())) {
          K = Kstar.children()[1];
        }
        else {
          continue;
        }
        // boost to B rest frame
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(B.momentum().betaVec());
        FourMomentum pChi   = boost.transform(chic.momentum());
        FourMomentum pKstar = boost.transform(Kstar.momentum());
        FourMomentum pK     = boost.transform(K    .momentum());
        // axes
        Vector3 axisX = pChi.p3().unit();
        // kaon helicity angle
        LorentzTransform boostK = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
        double cosK = -axisX.dot(boostK.transform(pK).p3().unit());
        // fill hists
        _h->fill(cosK);
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


  RIVET_DECLARE_PLUGIN(BELLE_2006_I689881);

}
