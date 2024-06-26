// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief kaon spetra at Upsilon(4S) and nearby continuum
  class CLEO_1982_I180612 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1982_I180612);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_c[ix],"TMP/c_"+toString(ix+1));
        book(_h_mult[ix], 1, 1, ix+1);
        for (unsigned int iy=0; iy<2; ++iy) {
          book(_h[ix][iy], 2, 1+iy, 1+ix);
        }
      }
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for (const Particle & p: mother.children()) {
        const int id = abs(p.pid());
        if (id == 310 || id == 130 || id == 321) {
          unstable.push_back(p);
        }
        else if(!p.children().empty()) {
          findDecayProducts(p, unstable);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==300553);
      // continuum
      if (upsilons.empty()) {
        _c[0]->fill();
        for (const Particle& p : ufs.particles(Cuts::pid==130 ||
                                               Cuts::pid==310 ||
                                               Cuts::abspid==321)) {
          int id = p.abspid();
          double modp = p.p3().mod();
          if (id==321) {
            _h[0][0]->fill(modp);
            _h_mult[0]->fill(_edges[0]);
          }
          else {
            _h[0][1]->fill(modp);
            _h_mult[0]->fill(_edges[1]);
          }
        }
      }
      /// found an upsilon
      else {
        for (const Particle& ups : upsilons) {
          _c[1]->fill();
          Particles unstable;
          // Find the decay products we want
          findDecayProducts(ups,unstable);
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 0.001) {
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          }
          for( const Particle & p : unstable) {
            int id = p.abspid();
            double modp = cms_boost.transform(p.momentum()).p3().mod();
            if (id==321) {
              _h[1][0]->fill(modp);
              _h_mult[1]->fill(_edges[0]);
            }
            else {
              _h[1][1]->fill(modp);
              _h_mult[1]->fill(_edges[1]);
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0;ix<2;++ix) {
        scale(_h_mult[ix], 1./ *_c[ix]);
        for (unsigned int iy=0;iy<2;++iy) {
          scale(_h[ix][iy],crossSection()/nanobarn/sumOfWeights());
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][2];
    BinnedHistoPtr<string> _h_mult[2];
    CounterPtr _c[2];
    vector<string> _edges = { "E+ E- --> (K+ + K-) X", "E+ E- --> K0 X" };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1982_I180612);

}
