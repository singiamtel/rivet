// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief ARGUS lepton spectra
  class ARGUS_1990_I298245 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1990_I298245);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // histograms
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix], 2,1,1+ix);
      }
      book(_c,"TMP/nb");
    }


    void findDecayProducts(const Particle& parent, Particles& ee, Particles& mu) {
      for (const Particle& p : parent.children()) {
      	if (p.abspid() == PID::EMINUS) {
          ee.push_back(p);
        }
        else if (p.abspid() == PID::MUON) {
          mu.push_back(p);
        }
      	else if (!p.children().empty()) {
      	  findDecayProducts(p,ee,mu);
      	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==300553)) {
      	// boost to rest frame
      	LorentzTransform cms_boost;
      	if (p.p3().mod() > 1*MeV) {
      	  cms_boost = LorentzTransform::mkFrameTransformFromBeta(p.mom().betaVec());
        }
      	// loop over decay products
      	for (const Particle& p2 : p.children()) {
      	  if (p2.abspid()==511 || p2.abspid()==521) {
      	    _c->fill();
            // find decay products
            Particles leptons[2];
            findDecayProducts(p2, leptons[0], leptons[1]);
            for (unsigned int ix=0; ix<2; ++ix) {
              for (const Particle& lep : leptons[ix]) {
                FourMomentum pl = cms_boost.transform(lep.mom());
                _h[ix]->fill(pl.p3().mod());
              }
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, 1.0/ *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1990_I298245);

}
