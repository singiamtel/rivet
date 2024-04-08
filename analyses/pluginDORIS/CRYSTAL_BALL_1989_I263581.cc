// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Electron spectrum  in B decays
  class CRYSTAL_BALL_1989_I263581 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CRYSTAL_BALL_1989_I263581);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // Book histograms
      // specify custom binning
      book(_h_all, 1, 1, 1);
      book(_nB, "/TMP/nB");
    }

    void findDecayProducts(Particle parent, Particles& em) {
      for (const Particle& p : parent.children()) {
        if(p.abspid() == PID::EMINUS) {
          em.push_back(p);
        }
        else  {
          findDecayProducts(p,em);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty())  _edges = _h_all->xEdges();
      // find and loop over Upslion(4S)
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==300553)) {
      	// boost to rest frame
      	LorentzTransform cms_boost;
      	if (p.p3().mod() > 1*MeV) {
      	  cms_boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
        }
      	// loop over decay products
      	for (const Particle& p2 : p.children()) {
      	  if (p2.abspid()==511 || p2.abspid()==521) {
            _nB->fill();
      	    Particles em;
            findDecayProducts(p2, em);
            for (const Particle& electron : em) {
              const double en = cms_boost.transform(electron.momentum()).E();
              _h_all->fill(map2string(en));
            }
      	  }
      	}
      }
    }

    string map2string(const double val) const {
      const size_t idx = _axis.index(val);
      if (idx && idx <= _edges.size())  return _edges[idx-1];
      return "OTHER";
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize to number of B decays
      scale(_h_all, 1./0.05/ *_nB);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_all;
    CounterPtr _nB;
    YODA::Axis<double> _axis = YODA::Axis<double>(48, 0.6, 3.0);
    vector<string> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CRYSTAL_BALL_1989_I263581);

}
