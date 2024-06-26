// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> p pbar K+, pi+
  class BELLE_2008_I754259 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2008_I754259);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      declare(BP, "BP");
      // histograms
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_mass [ix],1+ix,1,1);
        book(_h_angle[ix],3,1,1+ix);
      }
      book(_h_angle2, {2.*.9382722, 2.0, 2.2, 2.4, 2.6, 2.85});
      for (unsigned int ix=0; ix<5; ++ix) {
        book(_h_angle2->bin(ix+1), 4, 1, 1+ix);
      }
      book(_c,"TMP/nB");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for (unsigned int ix=0; ix<BP.decaying().size(); ++ix) {
        int imode=0,sign=1;
        _c->fill();
        if (BP.modeMatches(ix,3,mode1)) {
          imode= 0;
          sign = 1;
        }
        else if (BP.modeMatches(ix,3,mode1CC)) {
          imode= 0;
          sign =-1;
        }
        else if (BP.modeMatches(ix,3,mode2)) {
          imode= 1;
          sign = 1;
        }
        else if (BP.modeMatches(ix,3,mode2CC)) {
          imode= 1;
          sign =-1;
        }
        else {
          continue;
        }
        const Particle& pp   = BP.decayProducts()[ix].at( sign*2212)[0];
       	const Particle& pbar = BP.decayProducts()[ix].at(-sign*2212)[0];
        const Particle & meson = BP.decayProducts()[ix].at( sign*(imode==0 ? 321 : 211))[0];
        FourMomentum pbaryon = pp.mom()+pbar.mom();
        double mass = pbaryon.mass();
        _h_mass[imode]->fill(mass);
        if (mass>2.85) continue;
        // boost to B rest frame
        LorentzTransform boost =
          LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix]. mom().betaVec());
        FourMomentum pp1 = boost.transform(pbar.mom());
        pbaryon = boost.transform(pbaryon);
        Vector3 axis = boost.transform(meson.mom()).p3().unit();
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(pbaryon.betaVec());
        double cosp = axis.dot(boost1.transform(pp1).p3().unit());
        _h_angle[imode]->fill(cosp);
        if (imode==0) _h_angle2->fill(mass,cosp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_mass,  1e6/ *_c);
      scale(_h_angle, 1e6/ *_c);
      normalize(_h_angle2);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2],_h_angle[2];
    Histo1DGroupPtr _h_angle2;
    CounterPtr _c;
    const map<PdgId,unsigned int> mode1   = { { 2212,1}, {-2212,1}, { 321,1} };
    const map<PdgId,unsigned int> mode1CC = { { 2212,1}, {-2212,1}, {-321,1} };
    const map<PdgId,unsigned int> mode2   = { { 2212,1}, {-2212,1}, { 211,1} };
    const map<PdgId,unsigned int> mode2CC = { { 2212,1}, {-2212,1}, {-211,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2008_I754259);

}
