// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> e spectrum
  class CLEOII_2004_I647288 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_2004_I647288);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // histograms
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix],2,1,1+ix);
        book(_p[ix],1,1,1+ix);
      }
      book(_h[2],3,1,1);
      book(_c,"TMP/nb");
    }

    void findDecayProducts(const Particle& parent,
			   vector<Particles>&  em, vector<Particles>& ep,
			   vector<Particles>& nue, vector<Particles>& nueBar,
			   bool& charm, bool secondary) {
      for (const Particle& p : parent.children()) {
      	if (PID::isCharmHadron(p.pid())) {
      	  charm=true;
          findDecayProducts(p,em,ep,nue,nueBar,charm,true);
      	}
      	else if (p.pid() == PID::EMINUS || p.pid()==PID::MUON) {
      	  em[secondary].push_back(p);
      	}
      	else if (p.pid() == PID::EPLUS || p.pid()==PID::ANTIMUON) {
      	  ep[secondary].push_back(p);
      	}
      	else if (p.pid() == PID::NU_E  || p.pid()==PID::NU_MU) {
      	  nue[secondary].push_back(p);
      	}
      	else if (p.pid() == PID::NU_EBAR || p.pid()==PID::NU_MUBAR) {
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
      if(_edges.empty()) _edges = _p[0]->xEdges();
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
            _c->fill();
            // find decay products
            bool charm = false;
            vector<Particles> em(2),ep(2),nue(2),nueBar(2);
            findDecayProducts(p2,em,ep,nue,nueBar,charm,false);
            for (unsigned int isec=0;isec<2;++isec) {
              FourMomentum pl,pnu;
              if (em[isec].size()==1 && nueBar[isec].size()==1 && em[isec][0].pid()+1==-nueBar[isec][0].pid()) {
                pl  = cms_boost.transform(em[isec][0].momentum());
                pnu = cms_boost.transform(nueBar[isec][0].momentum());
              }
              else if (ep[isec].size()==1 && nue[isec].size()==1 && nue[isec][0].pid()==-ep[isec][0].pid()+1) {
                pl  = cms_boost.transform(ep[isec][0].momentum());
                pnu = cms_boost.transform(nue[isec][0].momentum());
              }
              else {
                continue;
              }
              _h[isec]->fill(pl.p3().mod());
              if (charm && isec==0) {
                _h[2]->fill(pl.p3().mod());
                double Emin=0.6;
                for(size_t i=0; i<_p[0]->numBins();++i) {
                  if(pl.E()>Emin) {
                    _p[0]->fill(_edges[i],pl.E());
                    _p[1]->fill(_edges[i],sqr(pl.E()));
                  }
                  Emin+=0.1;
                }
              }
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 0.5 from counting e and mu
      scale(_h, 0.5/ *_c);
      // compute <(el^2-<el>^2>
      BinnedEstimatePtr<string> tmp;
      book(tmp,1,1,3);
      for (const auto& b0 : _p[0]->bins()) {
        const auto& b1 = _p[1]->bin(b0.index());
        const double val = b1.yMean()-sqr(b0.yMean());
        const double err = val*sqrt(sqr(b1.relErrW())+4.*sqr(b0.relErrW()));
        tmp->bin(b0.index()).set(val, err);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    BinnedProfilePtr<string> _p[2];
    CounterPtr _c;
    vector<string> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_2004_I647288);

}
