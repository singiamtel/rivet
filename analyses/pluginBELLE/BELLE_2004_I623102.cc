// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> D** pi
  class BELLE_2004_I623102 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2004_I623102);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 411);
      BP.addStable(-411);
      BP.addStable( 413);
      BP.addStable(-413);
      declare(BP, "BP");
      // histograms
      book(_h_mass[0],1,1,1);
      book(_h_mass[1],4,1,1);
      book(_b_mass,{-1.,-0.67,-0.33,0,0.33,0.67,1.});
      for (unsigned int ix=0; ix<6; ++ix) {
        book(_b_mass->bin(ix+1), 2, 1, 1+ix);
      }
      for (unsigned int iy=0; iy<4; ++iy) {
        if (iy==0) book(_b_angle[0 ], {0.,5.,5.9,6.2,400});
        else       book(_b_angle[iy], {0.,5.76,5.98,6.15,400});
        for (unsigned int ix=0; ix<4; ++ix) {
          if (iy==0) book(_b_angle[0 ]->bin(ix+1),3, 1,1+ix);
          else       book(_b_angle[iy]->bin(ix+1),5,iy,1+ix);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      for (unsigned int ix=0; ix<BP.decaying().size(); ++ix) {
      	int sign = BP.decaying()[ix].pid()/BP.decaying()[ix].abspid();
        int imode=0;
      	if      (sign>0 && BP.modeMatches(ix,3,mode1))   imode=0;
      	else if (sign<0 && BP.modeMatches(ix,3,mode1CC)) imode=0;
      	else if (sign>0 && BP.modeMatches(ix,3,mode2))   imode=1;
      	else if (sign<0 && BP.modeMatches(ix,3,mode2CC)) imode=1;
        else continue;
        const Particles& pip = BP.decayProducts()[ix].at( sign*211);
        const Particle & Dm  = BP.decayProducts()[ix].at(-sign*(411+imode*2))[0];
        // find the mnimum Dpi mass
        double mDpi[2] = {(Dm.mom()+pip[0].mom()).mass(),
                          (Dm.mom()+pip[1].mom()).mass()};
        unsigned int iloc=0;
        if (mDpi[0]>mDpi[1]) {
          iloc=1;
          swap(mDpi[0],mDpi[1]);
        }
        _h_mass[imode]->fill(mDpi[0]);
        // compute the helicity angles
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix].mom().betaVec());
        // D pi mass and angle
        FourMomentum pDpi = Dm.mom()+pip[iloc].mom();
        pDpi = boost1.transform(pDpi);
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDpi.betaVec());
        FourMomentum pPi  = boost1.transform(pip[iloc].mom());
        FourMomentum pPi3 = boost2.transform(pPi);
        double cDpi = -pPi3.p3().unit().dot(pDpi.p3().unit());
        // fill histos
        if (imode==0) {
          _b_mass->fill(cDpi,mDpi[0]);
          _b_angle[0]->fill(sqr(mDpi[0]),cDpi);
        }
        else {
          _b_angle[1]->fill(sqr(mDpi[0]),cDpi);
          // find the pion from the D* decay
          Particle pi2;
          if (Dm.children().size()!=2) continue;
          if ((Dm.children()[0].pid()==-sign*411 ||
               Dm.children()[0].pid()==-sign*411) &&
              (Dm.children()[1].pid()==-sign*211 ||
               Dm.children()[1].pid()==      111)) {
            pi2 = Dm.children()[1];
          }
          else if ((Dm.children()[1].pid()==-sign*411 ||
                    Dm.children()[1].pid()==-sign*411) &&
                   (Dm.children()[0].pid()==-sign*211 ||
                    Dm.children()[0].pid()==      111)) {
            pi2 = Dm.children()[0];
          }
          else {
            continue;
          }
          FourMomentum pPi2 = boost1.transform(pi2.mom());
          FourMomentum pDm  = boost1.transform( Dm.mom());
          Vector3 axis = pDm.p3().unit();
          Vector3 trans1 = pPi .p3() - pPi .p3().dot(axis)*axis;
          Vector3 trans2 = pPi2.p3() - pPi2.p3().dot(axis)*axis;
          const double chi = atan2(trans1.cross(trans2).dot(axis),trans1.dot(trans2));
          _b_angle[3]->fill(sqr(mDpi[0]),chi);
          pDm  = boost2.transform(pDm );
          pPi2 = boost2.transform(pPi2);
          LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pDm.betaVec());
          pPi2 = boost3.transform(pPi2);
          const double cTheta = pPi2.p3().unit().dot(pPi3.p3().unit());
          _b_angle[2]->fill(sqr(mDpi[0]),cTheta);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_mass, 1.0, false);
      normalize(_b_mass, 1.0, false);
      normalize(_b_angle, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2];
    Histo1DGroupPtr _b_mass, _b_angle[4];
    const map<PdgId,unsigned int> mode1   = { { 211,2}, {-411,1}};
    const map<PdgId,unsigned int> mode1CC = { {-211,2}, { 411,1}};
    const map<PdgId,unsigned int> mode2   = { { 211,2}, {-413,1}};
    const map<PdgId,unsigned int> mode2CC = { {-211,2}, { 413,1}};
    /// @}
  };


  RIVET_DECLARE_PLUGIN(BELLE_2004_I623102);

}
