// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  class LHCF_2012_I1115479 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(LHCF_2012_I1115479);

  public:

    void init() {
      declare(UnstableParticles(),"UFS");

      book(_g, {8.9, 9., 9.2, 9.4, 9.6, 10., 11.});
      for (auto& b : _g->bins()) {
        book(b, b.index(), 1, 1);
      }
    }


    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      const double dphi = TWOPI;

      for (const Particle& p : ufs.particles()) {
        if (p.pid() == 111) {
          double pT = p.pT();
          double y  = p.rapidity();

          if (pT > 0.6*GeV) continue;

          const double scaled_weight = 1.0/(dphi*pT/GeV);
          _g->fill(y, pT/GeV, scaled_weight);
        }
      }
    }


    void finalize() {
      scale(_g, 1./sumOfWeights());
      divByGroupWidth(_g);
    }

  private:

    Histo1DGroupPtr _g;

  };


  RIVET_DECLARE_PLUGIN(LHCF_2012_I1115479);

}
