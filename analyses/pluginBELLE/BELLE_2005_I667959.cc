// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> K pi pi K K K
  class BELLE_2005_I667959 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2005_I667959);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable(PID::PI0);
      BP.addStable(PID::K0S);
      BP.addStable(PID::ETA);
      BP.addStable(PID::ETAPRIME);
      declare(BP, "BP");
      vector<double> bins = {0.,5.,10.,15.,20.,30.};
      vector<double> bins2 = {0.,1.05,5.,10.,30.};
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_mass3[ix], 2, 1, 1+ix);
        book(_h_mass2[0][ix],bins);
        if (ix==0) book(_h_mass2[1][ix],bins );
        else       book(_h_mass2[1][ix],bins2);
        for(unsigned int iy=0; iy<3; ++iy) {
          book(_h_mass[ix][iy],1+4*ix,1,1+iy);
          if (iy==2) continue;
          book(_h_mass_zoom[ix][iy],1+4*ix,2,1+iy);
        }
        for(unsigned int iz=0;iz<5;++iz) {
          book(_h_mass2[0][ix]->bin(1+iz),3,1+ix,1+iz);
          if (ix==1 &&iz==4) continue;
          book(_h_mass2[1][ix]->bin(iz+1),6,1+ix,1+iz);
        }
      }
      for (unsigned int ix=0;ix<3;++ix) {
        vector<double> bins3;
        if (ix==0)      bins3 = {0.82,0.97,1.,1.76};
        else if (ix==1) bins3 = {0.,0.9,1.06};
        else if (ix==2) bins3 = {0.,1.05,3.0};
        book(_h_angle[ix],bins3);
        for (unsigned int iy=0; iy<2; ++iy) {
          if      (ix==0) book(_h_angle[ix]->bin(2*iy+1),4,1,1+iy);
          else if (ix==1) book(_h_angle[ix]->bin(iy+1)  ,4,1,3+iy);
          else if (ix==2) book(_h_angle[ix]->bin(iy+1)  ,7,1,1+iy);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const double mD0 = 1.86484, mJpsi = 3.096916, mPsi2S=3.6861;
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for (unsigned int ix=0; ix<BP.decaying().size(); ++ix) {
        int sign=BP.decaying()[ix].pid()/BP.decaying()[ix].abspid();
      	// boost to B rest frame
        LorentzTransform boost =
                LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix]. mom().betaVec());
        if (BP.modeMatches(ix,3,mode1) || BP.modeMatches(ix,3,mode1CC)) {
          const Particle & Kp   = BP.decayProducts()[ix].at( sign*321)[0];
          const Particle & pip  = BP.decayProducts()[ix].at( sign*211)[0];
          const Particle & pim  = BP.decayProducts()[ix].at(-sign*211)[0];
          FourMomentum pKpi = Kp.mom() +pim.mom();
          double mKpi  = pKpi.mass();
          if (abs(mKpi -mD0)<0.1 ) continue;
          FourMomentum ppipi = pip.mom()+pim.mom();
          double mpipi = ppipi.mass();
          if (abs(mpipi -mD0)<0.015 || abs(mpipi-mJpsi)<0.07 || abs(mpipi-mPsi2S)<0.05) continue;
          _h_mass[0][0]->fill(mKpi);
          _h_mass[0][1]->fill(mpipi);
          _h_mass[0][2]->fill((Kp.mom()+pip.mom()).mass());
          _h_mass_zoom[0][0]->fill(mKpi);
          _h_mass_zoom[0][1]->fill(mpipi);
          if (sqr(mpipi)<11.) _h_mass3[0]->fill(mKpi);
          else                _h_mass3[1]->fill(mKpi);
          _h_mass2[0][0]->fill(sqr(mpipi),mKpi );
          _h_mass2[0][1]->fill(sqr(mKpi ),mpipi);
          // helicity angles
          // Kpi
          pKpi = boost.transform(pKpi);
          LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKpi.betaVec());
          FourMomentum ppi = boost2.transform(boost.transform(pim.mom()));
          const double cK = ppi.p3().unit().dot(pKpi.p3().unit());
          _h_angle[0]->fill(mKpi,cK);
          // pi pi
          ppipi = boost.transform(ppipi);
          LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(ppipi.betaVec());
          ppi = boost3.transform(boost.transform(pim.mom()));
          const double cPi = ppi.p3().unit().dot(ppipi.p3().unit());
          _h_angle[1]->fill(mpipi,cPi);
        }
        else if (BP.modeMatches(ix,3,mode2) || BP.modeMatches(ix,3,mode2CC)) {
          const Particles& Kp   = BP.decayProducts()[ix].at( sign*321);
          const Particle & Km   = BP.decayProducts()[ix].at(-sign*321)[0];
          FourMomentum pKK[2] = {Kp[0].mom()+Km.mom(), Kp[1].mom()+Km.mom()};
          double mKK[2] = {pKK[0].mass(),pKK[1].mass()};
          if (mKK[0]>mKK[1]) {
            swap(mKK[0],mKK[1]);
            swap(pKK[0],pKK[1]);
          }
          _h_mass[1][0]->fill(mKK[0]);
          _h_mass[1][1]->fill(mKK[1]);
          _h_mass[1][2]->fill((Kp[0].mom()+Kp[1].mom()).mass());
          _h_mass_zoom[1][0]->fill(mKK[0]);
          _h_mass_zoom[1][1]->fill(mKK[1]);
          _h_mass2[1][0]->fill(sqr(mKK[1]),mKK[0]);
          _h_mass2[1][1]->fill(sqr(mKK[0]),mKK[1]);
          // helicity angles
          // Kpi
          pKK[0] = boost.transform(pKK[0]);
          LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKK[0].betaVec());
          FourMomentum pK = boost2.transform(boost.transform(Km.mom()));
          double cK = pK.p3().unit().dot(pKK[0].p3().unit());
          _h_angle[2]->fill(mKK[0],cK);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0; ix<2; ++ix) {
        normalize(_h_mass3[ix], 1.0, false);
        normalize(_h_mass[ix], 1.0, false);
        normalize(_h_mass_zoom[ix], 1.0, false);
        normalize(_h_mass2[ix], 1.0, false);
      }
      normalize(_h_angle, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2][3],_h_mass_zoom[2][2],_h_mass3[2];
    Histo1DGroupPtr _h_mass2[2][2],_h_angle[3];
    const map<PdgId,unsigned int> mode1   = { { 321,1}, { 211,1}, {-211,1} };
    const map<PdgId,unsigned int> mode1CC = { {-321,1}, { 211,1}, {-211,1} };
    const map<PdgId,unsigned int> mode2   = { { 321,2}, {-321,1} };
    const map<PdgId,unsigned int> mode2CC = { {-321,2}, { 321,1} };
  /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2005_I667959);

}
