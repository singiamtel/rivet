// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {

  class CMS_2015_I1342266 : public Analysis {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2015_I1342266);

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h_dy_pT[ix], {0.0,0.6,1.2});
        for (unsigned int iy=1; iy<3; ++iy) {
          book(_h_dy_pT[ix]->bin(iy), iy, 1, 1+ix);
        }
        book(_h_pT[ix], 3, 1, 1+ix);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==553 || Cuts::pid==100553 || Cuts::pid==200553)) {
        const double absrap = p.absrap();
        const double xp = p.perp();
        if (absrap>1.2) continue;
        unsigned int iups = p.pid()/100000;
        _h_pT   [iups]->fill(xp);
        _h_dy_pT[iups]->fill(absrap,xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double factor = crossSection() / femtobarn/ sumOfWeights();
      const vector<double> brs = {0.0248,0.0193,0.0218};
      for (unsigned int ix=0; ix<3; ++ix) {
        scale(_h_dy_pT[ix], factor*brs[ix]); //branching ratio
        divByGroupWidth(_h_dy_pT[ix]);
        scale(_h_pT[ix],factor * brs[ix]);
      }
      for (unsigned int i=1; i<3; ++i) {
        Estimate1DPtr tmp;
        book(tmp,i+3,1,1);
        divide(_h_dy_pT[1]->bin(i), _h_dy_pT[0]->bin(i), tmp);
        book(tmp,i+3,1,2);
        divide(_h_dy_pT[2]->bin(i), _h_dy_pT[0]->bin(i), tmp);
      }
      Estimate1DPtr tmp;
      book(tmp,6,1,1);
      divide(_h_pT[1], _h_pT[0], tmp);
      book(tmp,6,1,2);
      divide(_h_pT[2], _h_pT[0], tmp);
    }

  private:

    Histo1DGroupPtr _h_dy_pT[3];
    Histo1DPtr _h_pT[3];

  };

  RIVET_DECLARE_PLUGIN(CMS_2015_I1342266);
}
