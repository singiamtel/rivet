// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief BS0 -> J/psi f0
  class LHCB_2013_I1081268 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2013_I1081268);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==531);
      declare(ufs, "UFS");
      DecayedParticles BS0(ufs);
      BS0.addStable( 443);
      declare(BS0, "BS0");
      // histograms
      for (unsigned int ix=0;ix<2;++ix) {
        book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BS0 = apply<DecayedParticles>(event, "BS0");
      static const map<PdgId,unsigned int> mode = { { 211,1}, {-211,1}, { 443,1}};
      for (unsigned int ix=0; ix<BS0.decaying().size(); ++ix) {
       	if (!BS0.modeMatches(ix,3,mode)) continue;
      	const Particle& pip  = BS0.decayProducts()[ix].at( 211)[0];
      	const Particle& pim  = BS0.decayProducts()[ix].at(-211)[0];
      	const Particle& JPsi = BS0.decayProducts()[ix].at( 443)[0];
        FourMomentum ppipi = pip.mom()+pim.mom();
        // mass cut and presence of f_0(980)
        const double mpipi = ppipi.mass();
        if(abs(mpipi-0.98)>0.9) continue;
        if(count(BS0.decaying()[ix].children(), hasPID(9010221))!=1) continue;
        // f_0 helicity angle
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BS0.decaying()[ix].mom().betaVec());
        FourMomentum ppsi  = boost1.transform(JPsi.mom());
        Vector3 axis1 = ppsi.p3().unit();
        ppipi = boost1.transform(ppipi);
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(ppipi.betaVec());
        FourMomentum ppip = boost2.transform(boost1.transform(pip.mom()));
        _h[1]->fill(ppip.p3().unit().dot(axis1));
        // helicity angle find J.psi leptonic children
        if (JPsi.children().size()!=2) continue;
        if (JPsi.children()[0].pid()!=-JPsi.children()[1].pid()) continue;
        if (JPsi.children()[0].abspid()!=PID::EMINUS &&
            JPsi.children()[0].abspid()!=PID::MUON) continue;
        Particle lm = JPsi.children()[0];
        Particle lp = JPsi.children()[1];
        if (lm.pid()<0) swap(lm,lp);
        LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(ppsi.betaVec());
        FourMomentum plp = boost3.transform(boost1.transform(lp.mom()));
        _h[0]->fill(-plp.p3().unit().dot(axis1));
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
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2013_I1081268);

}
