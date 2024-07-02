// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> e spectrum
  class ARGUS_1993_I357133 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1993_I357133);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // histograms
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix],1+ix,1,1);
      }
      book(_c,"TMP/nb");
    }

    void findDecayProducts(const Particle& parent,
                           vector<Particles>&  em, vector<Particles>& ep,
                           vector<Particles>& nue, vector<Particles>& nueBar,
                           bool& charm, bool secondary) {
      for (const Particle & p : parent.children()) {
      	if (PID::isCharmHadron(p.pid())) {
      	  charm=true;
          findDecayProducts(p,em,ep,nue,nueBar,charm,true);
      	}
      	else if (p.pid() == PID::EMINUS) {
      	  em[secondary].push_back(p);
      	}
      	else if (p.pid() == PID::EPLUS) {
      	  ep[secondary].push_back(p);
      	}
      	else if (p.pid() == PID::NU_E ) {
      	  nue[secondary].push_back(p);
      	}
      	else if (p.pid() == PID::NU_EBAR ) {
      	  nueBar[secondary].push_back(p);
      	}
      	else if (PID::isBottomHadron(p.pid())) {
      	  findDecayProducts(p,em,ep,nue,nueBar,charm,false);
      	}
      	else if (!PID::isHadron(p.pid())) {
      	  findDecayProducts(p,em,ep,nue,nueBar,charm,secondary);
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
      	for (const Particle & p2 : p.children()) {
      	  if (p2.abspid()==511 || p2.abspid()==521) {
      	    _c->fill();
            // find decay products
            bool charm = false;
            vector<Particles> em(2),ep(2),nue(2),nueBar(2);
            findDecayProducts(p2,em,ep,nue,nueBar,charm,false);
            for (unsigned int isec=0; isec<2; ++isec) {
              FourMomentum pl,pnu;
              if (em[isec].size()==1 && nueBar[isec].size()==1 && em[isec][0].pid()+1==-nueBar[isec][0].pid()) {
                pl  = cms_boost.transform(em[isec][0].mom());
                pnu = cms_boost.transform(nueBar[isec][0].mom());
              }
              else if (ep[isec].size()==1 && nue[isec].size()==1 && nue[isec][0].pid()==-ep[isec][0].pid()+1) {
                pl  = cms_boost.transform(ep[isec][0].mom());
                pnu = cms_boost.transform(nue[isec][0].mom());
              }
              else {
                continue;
              }
              _h[isec]->fill(pl.p3().mod());
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


  RIVET_DECLARE_PLUGIN(ARGUS_1993_I357133);

}
