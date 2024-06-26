// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e- spectrum in semileptonic B decays
  class BELLE_2005_I665011 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2005_I665011);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(Cuts::abspid==511 ||
				Cuts::abspid==521),"UFS");
      // histos
      for (unsigned int ix=0;ix<2;++ix) {
        book(_h[ix],1,1,1+ix);
        book(_c[ix],"TMP/c_"+toString(ix+1));
      }
    }

    void findDecayProducts(const Particle& parent, Particles& em, Particles& ep,
                           Particles& nue, Particles& nueBar, bool& charm) {
      for (const Particle& p : parent.children()) {
        if (PID::isCharmHadron(p.pid())) {
          charm=true;
        }
        else if (p.pid() == PID::EMINUS) {
          em.push_back(p);
        }
        else if (p.pid() == PID::EPLUS) {
          ep.push_back(p);
        }
        else if (p.pid() == PID::NU_E) {
          nue.push_back(p);
        }
        else if (p.pid() == PID::NU_EBAR) {
          nueBar.push_back(p);
        }
        else if (PID::isBottomHadron(p.pid())) {
          findDecayProducts(p,em,ep,nue,nueBar,charm);
        }
        else if (!PID::isHadron(p.pid())) {
          findDecayProducts(p,em,ep,nue,nueBar,charm);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find and loop over B mesos
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        unsigned int iB = (p.abspid()-511)/10;
        _c[iB]->fill();
        // find decay products
        bool charm = false;
        Particles em,ep,nue,nueBar;
        findDecayProducts(p,em,ep,nue,nueBar,charm);
        if (!charm) continue;
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
        double pmod(0.);
        if (em.size()==1 && nueBar.size()==1) {
          pmod = boost.transform(em[0].momentum()).p3().mod();
        }
        else if (ep.size()==1 && nue.size()==1) {
          pmod = boost.transform(ep[0].momentum()).p3().mod();
        }
        else {
          continue;
        }
        _h[iB]->fill(pmod);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize hists
      for (unsigned int ix=0; ix<2; ++ix) {
	      scale(_h[ix],1./ *_c[ix]);
      }
      // ratio
      Estimate1DPtr ratio;
      book(ratio,1,1,3);
      divide(_h[1], _h[0], ratio);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    CounterPtr _c[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2005_I665011);

}
