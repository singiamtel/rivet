// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi spectrum in Upsilon(1S) decay
  class CLEOIII_2004_I654639 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOIII_2004_I654639);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==553), "UFS");
      // histos
      for (unsigned int ix=0;ix<2;++ix) {
        book(_h[ix],1,1,1+ix);
       }
    }

   /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for (const Particle& p: mother.children()) {
        const int id = abs(p.pid());
        if (id == 443) {
          unstable.push_back(p);
        }
        if (!p.children().empty()) {
          findDecayProducts(p, unstable);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      for (const Particle& ups : apply<UnstableParticles>(event, "UFS").particles()) {
        // boost to rest frame
        LorentzTransform boost;
        if (ups.p3().mod() > 1*MeV) {
          boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
        }
        // find  psi
        Particles unstable;
        findDecayProducts(ups,unstable);
        // loop over psi/psi(2S)
        for (const Particle& psi : unstable) {
          double Pmax = 0.5*(sqr(ups.mass())-sqr(psi.mass()))/ups.mass();
          const FourMomentum p2 = boost.transform(psi.momentum());
          const double xP = p2.p3().mod()/Pmax;
          _h[0]->fill(xP);
          _h[1]->fill(xP);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, crossSection()/sumOfWeights()/picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOIII_2004_I654639);

}
