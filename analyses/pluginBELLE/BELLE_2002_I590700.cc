// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> D(*) K K(*)
  class BELLE_2002_I590700 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2002_I590700);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 or Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 411);
      BB.addStable(-411);
      BB.addStable( 421);
      BB.addStable(-421);
      BB.addStable( 413);
      BB.addStable(-413);
      BB.addStable( 423);
      BB.addStable(-423);
      BB.addStable( 310);
      BB.addStable( 313);
      BB.addStable(-313);
      declare(BB, "BB");
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_mass[ix],1+2*ix,1,1);
        for (unsigned int iy=0;iy<2;++iy) {
          book(_h_angle[ix][iy],2,1+iy,1+ix);
        }
      }
      book(_h_angle[2][0],3,2,1);
    }

    bool isK(int pid) const {
      return pid==130 || pid==310 || pid==311 || pid==321;
    }

    bool isPi(int pid) const {
      return pid==211 || pid==111;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const double mDs = 1.96835;
      static const map<PdgId,unsigned int> modes  [5] = { { {-421,1}, { 321,1}, {-313,1} },
                                                          { {-411,1}, { 321,1}, {-313,1} },
                                                          { {-423,1}, { 321,1}, {-313,1} },
                                                          { {-413,1}, { 321,1}, {-313,1} },
                                                          { {-421,1}, { 321,1}, { 310,1} }};
      static const map<PdgId,unsigned int> modesCC[5] = { { { 421,1}, {-321,1}, { 313,1} },
                                                          { { 411,1}, {-321,1}, { 313,1} },
                                                          { { 423,1}, {-321,1}, { 313,1} },
                                                          { { 413,1}, {-321,1}, { 313,1} },
                                                          { { 421,1}, {-321,1}, { 310,1} }};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      for (unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
        int imode=-1;
        for (unsigned int iy=0; iy<5; ++iy) {
          if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,modes[iy])) {
            imode=iy;
            break;
          }
          else if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,modesCC[iy])) {
            imode=iy;
            break;;
          }
        }
        if (imode<0) continue;
        int sign = BB.decaying()[ix].pid()/BB.decaying()[ix].abspid();
        LorentzTransform boostB = LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix].mom().betaVec());
        if (imode<4) {
          const Particle& Km = BB.decayProducts()[ix].at( sign*321)[0];
          const Particle& K0 = BB.decayProducts()[ix].at(-sign*313)[0];
          FourMomentum pKK = boostB.transform(Km.mom()+K0.mom());
          const double mKK=pKK.mass();
          if (abs(mKK-mDs)<0.02) continue;
          _h_mass[0]->fill(mKK);
          LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKK.betaVec());
          FourMomentum pK0 = boost2.transform(boostB.transform(K0.mom()));
          _h_angle[imode/2][0]->fill(pK0.p3().unit().dot(pKK.p3().unit()));
          // find Kstar decay products
          Particle KK;
          if (isK (K0.children()[0].abspid()) &&
              isPi(K0.children()[1].abspid())) {
            KK = K0.children()[0];
          }
          else if (isK(K0.children()[1].abspid()) &&
                  isPi(K0.children()[0].abspid())) {
            KK = K0.children()[1];
          }
          else {
            continue;
          }
          LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pK0.betaVec());
          FourMomentum pK = boost3.transform(boost2.transform(boostB.transform(KK.mom())));
          _h_angle[imode/2][1]->fill(pKK.p3().unit().dot(pK0.p3().unit()));
        }
        else {
          const Particle& Km = BB.decayProducts()[ix].at( sign*321)[0];
          const Particle& K0 = BB.decayProducts()[ix].at(      310)[0];
          FourMomentum pKK = boostB.transform(Km.mom()+K0.mom());
          const double mKK=pKK.mass();
          if (abs(mKK-mDs)<0.02) continue;
          _h_mass[1]->fill(mKK);
          LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKK.betaVec());
          FourMomentum pK0 = boost2.transform(boostB.transform(K0.mom()));
          _h_angle[2][0]->fill(pK0.p3().unit().dot(pKK.p3().unit()));
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0; ix<2; ++ix) {
        normalize(_h_mass[ix], 1.0, false);
        for (unsigned int iy=0; iy<2; ++iy) {
          normalize(_h_angle[ix][iy], 1.0, false);
        }
      }
      normalize(_h_angle[2][0],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2],_h_angle[3][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2002_I590700);

}
