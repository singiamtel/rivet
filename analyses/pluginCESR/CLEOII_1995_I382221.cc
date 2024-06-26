// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  J/Psi, Psi(2S), chi_c1 spectra at the Upsilon(4S)
  class CLEOII_1995_I382221 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1995_I382221);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::pid==300553), "UFS");
      // histos
      for (unsigned int ix=0; ix<4; ++ix) {
        book(_h[ix], 1+ix,1,1);
      }
      book(_c,"TMP/nUps");
    }

    void findDecayProducts(Particle parent, Particles & primary,
			   Particles & secondary, bool second) {
      for (const Particle& p :parent.children()) {
        if (p.pid()==443) {
          if (second) secondary.push_back(p);
          else        primary.push_back(p);
       	  continue;
       	}
       	else if (p.pid()==100443 || p.pid()==20443) {
          if (second) secondary.push_back(p);
          else        primary.push_back(p);
          findDecayProducts(p,primary,secondary,true);
        }
        else if (p.pid()==10441||p.pid()==445) {
          findDecayProducts(p,primary,secondary,true);
        }
        else if(!p.children().empty()) {
          findDecayProducts(p,primary,secondary,second);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        _c->fill();
        const LorentzTransform boost =
              LorentzTransform::mkFrameTransformFromBeta(p.mom().betaVec());
       	Particles primary,secondary;
       	findDecayProducts(p,primary,secondary,false);
       	for (const Particle& pp : primary) {
      	  const double pcm = boost.transform(pp.mom()).p3().mod();
          if (pp.pid()==443) {
            _h[0]->fill(pcm);
            _h[1]->fill(pcm);
          }
          else if(pp.pid()==100443) _h[2]->fill(pcm);
          else if(pp.pid()==20443)  _h[3]->fill(pcm);
        }
       	for (const Particle& pp : secondary) {
      	  const double pcm = boost.transform(pp.mom()).p3().mod();
          if (pp.pid()==443)        _h[0]->fill(pcm);
          else if (pp.pid()==100443) _h[2]->fill(pcm);
          else if (pp.pid()==20443)  _h[3]->fill(pcm);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 100 convert to % and 0.5 as 2 B decays per Upslion(4S)
	    scale(_h, 50./ *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1995_I382221);

}
