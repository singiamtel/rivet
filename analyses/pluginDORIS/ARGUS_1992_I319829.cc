// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s+ in B or continuum
  class ARGUS_1992_I319829 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1992_I319829);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      declare(Beam(), "Beams");
      // histos
      book(_n_Ds, 1, 1, 1);
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_Ds[ix], 2+ix,1,1);
      }
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(const Particle& mother, Particles& unstable) {
      for(const Particle& p: mother.children()) {
        const int id = p.abspid();
        if (id == 431) {
          unstable.push_back(p);
        }
        else if (!p.children().empty()) {
          findDecayProducts(p, unstable);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = 0.5*(beams.first.p3().mod() + beams.second.p3().mod());
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // Find the Upsilon(4S) among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==300553);
      // Continuum
      if (upsilons.empty()) {
        for (const Particle &p : ufs.particles(Cuts::abspid==431)) {
          const double xp = p.p3().mod()/sqrt(sqr(meanBeamMom)-sqr(p.mass()));
          _h_Ds[0]->fill(xp);
          _n_Ds->fill(10);
        }
      }
      else {
        for (const Particle& ups : upsilons) {
          Particles unstable;
          // Find the decay products we want
          findDecayProducts(ups, unstable);
          // boost to rest frame (if required)
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 1*MeV) {
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.mom().betaVec());
          }
          for(const Particle& p : unstable) {
            const FourMomentum p2 = cms_boost.transform(p.mom());
            double xp = p2.p3().mod()/sqrt(0.25*sqr(ups.mass())-sqr(p.mass()));
            _h_Ds[1]->fill(xp);
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // bre d+s -> phi pi+
      const double br = 0.045;
      normalize(_h_Ds, 1.0, false);
      scale(_n_Ds, br/sumOfWeights()*crossSection()/picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _n_Ds;
    Histo1DPtr _h_Ds[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1992_I319829);

}
