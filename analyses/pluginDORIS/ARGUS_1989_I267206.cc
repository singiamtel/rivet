// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Bbar0 -> D* l nu
  class ARGUS_1989_I267206 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1989_I267206);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(411); BB.addStable(-411);
      BB.addStable(421); BB.addStable(-421);
      BB.addStable(413); BB.addStable(-413);
      BB.addStable(423); BB.addStable(-423);
      BB.addStable(PID::PI0);
      declare(BB, "BB");
      // histos
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h[ix], 1+ix, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over B mesons
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for (unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
      	for (unsigned int il=0; il<2; ++il) {
          if (BB.modeMatches(ix,3,mode[il])) {
            int sign = il%2==0 ? -1 : 1;
            const Particle & DD = BB.decayProducts()[ix].at(sign*413)[0];
            const Particle & em = BB.decayProducts()[ix].at(sign*11)[0];
            FourMomentum qq  = BB.decaying()[ix].mom() - DD.mom();
            _h[1]->fill(em.mom().p3().mod());
            _h[2]->fill(qq.mass2());
            if(DD.children().size()!=2) continue;
            // boost to B rest frame
            LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix].mom().betaVec());
            FourMomentum pDD = boost1.transform(DD.mom());
            Particle D0;
            if (DD.children()[0].pid() ==sign*211 && DD.children()[1].pid()==sign*421) {
              D0 = DD.children()[1];
            }
            else if (DD.children()[1].pid()==sign*211 && DD.children()[0].pid()==sign*421) {
              D0 = DD.children()[0];
            }
            else {
              continue;
            }
            FourMomentum pD0 = boost1.transform(D0.mom());
            LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pDD.betaVec());
            pD0 = boost3.transform(pD0);
            double cTheta = pD0.p3().unit().dot(pDD.p3().unit());
            _h[0]->fill(cTheta);
          }
        }
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
    const map<PdgId,unsigned int> mode[2] = {{ {-413,1}, {-11,1}, { 12,1}}, { { 413,1}, { 11,1}, {-12,1}}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1989_I267206);

}
