// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> pi+ pi0 pi0
  class BELLE_2005_I651511 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2005_I651511);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 111);
      declare(BP, "BP");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix], 1+ix, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      for (unsigned int ix=0;ix<BP.decaying().size();++ix) {
      	int sign = 1;
      	if (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,3,mode))        sign= 1;
      	else if (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,3,modeCC)) sign=-1;
      	else continue;
       	const Particle& pip = BP.decayProducts()[ix].at( sign*211)[0];
       	const Particles& pi0 = BP.decayProducts()[ix].at(      111);
        FourMomentum pRho[2] = {pip.mom()+pi0[0].mom(),	pip.mom()+pi0[1].mom()};
        double mRho[2] = {pRho[0].mass(),pRho[1].mass()};
        unsigned int iloc = mRho[0]<mRho[1] ? 0 : 1;
        _h[1]->fill(mRho[iloc]);
        // compute the helicity angles
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix].mom().betaVec());
        FourMomentum prho = boost1.transform(pRho[iloc]);
        FourMomentum ppi  = boost1.transform(pip.mom());
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(prho.betaVec());
        ppi = boost2.transform(ppi);
        const double cTheta = ppi.p3().unit().dot(prho.p3().unit());
        _h[0]->fill(cTheta);
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
    const map<PdgId,unsigned int> mode   = { { 211,1}, { 111,2}};
    const map<PdgId,unsigned int> modeCC = { {-211,1}, { 111,2}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2005_I651511);

}
