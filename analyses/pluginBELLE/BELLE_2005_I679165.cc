// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> p pbar K or p Lambdabar pi
  class BELLE_2005_I679165 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2005_I679165);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(310);
      BB.addStable( 3122);
      BB.addStable(-3122);
      declare(BB, "BB");
      // histograms
      for (unsigned int ix=0; ix<4; ++ix) {
        book(_h_angle[ix], 2, 1, 1+ix);
        if (ix==3) continue;
        book(_h_mass[ix], 1, 1, 1+ix);
      }
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_c[ix],"TMP/nB_"+toString(ix+1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for (unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
        int sign = BB.decaying()[ix].pid()/BB.decaying()[ix].abspid();
        int imode=-1;
        if (BB.decaying()[ix].abspid()==511) _c[0]->fill();
        else                                 _c[1]->fill();
        if ((sign>0 &&BB.modeMatches(ix,3,mode1)) ||
            (sign<0 &&BB.modeMatches(ix,3,mode1CC))) imode = 0;
        else if ((sign>0 &&BB.modeMatches(ix,3,mode2)) ||
                 (sign<0 &&BB.modeMatches(ix,3,mode2CC))) imode = 1;
        else if ((sign>0 &&BB.modeMatches(ix,3,mode3)) ||
                 (sign<0 &&BB.modeMatches(ix,3,mode3CC))) imode = 2;
        else {
          continue;
        }
        // boost to B rest frame
        LorentzTransform boost =
          LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix]. mom().betaVec());
        if (imode<2) {
          const Particle& pp    = BB.decayProducts()[ix].at( sign*2212)[0];
          const Particle& pbar  = BB.decayProducts()[ix].at(-sign*2212)[0];
          const Particle& meson = BB.decayProducts()[ix].at( imode==0 ? sign*321 : 310)[0];
          FourMomentum pbaryon = pp.mom()+pbar.mom();
          double mass = pbaryon.mass();
          if (2.850<mass && mass<3.12 ) continue;
          if (3.315<mass && mass<3.735) continue;
          _h_mass[imode]->fill(mass);
          if (mass>2.85) continue;
          FourMomentum pp1 = boost.transform(pp.mom());
          pbaryon = boost.transform(pbaryon);
          Vector3 axis = boost.transform(meson.mom()).p3().unit();
          LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(pbaryon.betaVec());
          double cosp = axis.dot(boost1.transform(pp1).p3().unit());
          _h_angle[imode]->fill(cosp);
        }
        else {
          const Particle& pp    = BB.decayProducts()[ix].at( sign*2212)[0];
          const Particle& lbar  = BB.decayProducts()[ix].at(-sign*3122)[0];
          const Particle& meson = BB.decayProducts()[ix].at(-sign* 211)[0];
          FourMomentum pbaryon = pp.mom()+lbar.mom();
          double mass = pbaryon.mass();
          _h_mass[imode]->fill(mass);
          FourMomentum pp1 = boost.transform(pp.mom());
          FourMomentum pppi = pp.mom()+meson.mom();
          LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pppi.betaVec());
          Vector3 axis = boost.transform(lbar.mom()).p3().unit();
          double cosp = axis.dot(boost2.transform(pp1).p3().unit());
          _h_angle[3]->fill(cosp);
          if (mass>2.85) continue;
          pbaryon = boost.transform(pbaryon);
          axis = boost.transform(meson.mom()).p3().unit();
          LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(pbaryon.betaVec());
          cosp = axis.dot(boost1.transform(pp1).p3().unit());
          _h_angle[imode]->fill(cosp);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_mass [0], 1e6/ *_c[1]);
      scale(_h_mass [1], 1e6/ *_c[0]);
      scale(_h_mass [2], 1e6/ *_c[0]);
      scale(_h_angle[0], 1e6/ *_c[1]);
      for (unsigned int ix=1; ix<4; ++ix) {
        scale(_h_angle[ix], 1e6/ *_c[0]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[3], _h_angle[4];
    CounterPtr _c[2];
    const map<PdgId,unsigned int> mode1   = { { 2212,1}, {-2212,1}, { 321,1} };
    const map<PdgId,unsigned int> mode1CC = { { 2212,1}, {-2212,1}, {-321,1} };
    const map<PdgId,unsigned int> mode2   = { { 2212,1}, {-2212,1}, { 310,1} };
    const map<PdgId,unsigned int> mode2CC = { { 2212,1}, {-2212,1}, { 310,1} };
    const map<PdgId,unsigned int> mode3   = { { 2212,1}, {-3122,1}, {-211,1} };
    const map<PdgId,unsigned int> mode3CC = { {-2212,1}, { 3122,1}, { 211,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2005_I679165);

}
