// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B-> baryons
  class ARGUS_1992_I334962 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1992_I334962);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // histos
      for (unsigned int ix=0; ix < 2; ++ix) {
        book(_h_spect[ix],1+ix,1,1);
        book(_h_flav [ix],3,1,1+ix);
      }
    }

    void findDecayProducts(const Particle& parent, Particles& protons_all, Particles& protons_charm,
                           Particles& protons_had, Particles & lambda, bool cBaryon=false) {
      for (const Particle& p : parent.children()) {
       	if (p.abspid() == PID::PROTON) {
       	  protons_all.push_back(p);
          if (cBaryon) protons_charm.push_back(p);
          else         protons_had  .push_back(p);
        }
        else if (p.abspid()== PID::LAMBDA) {
          lambda.push_back(p);
          continue;
        }
        else if (!p.children().empty()) {
          bool cBaryonDec = cBaryon || ( PID::isCharmHadron(p.pid()) && p.abspid()%2==0);
          findDecayProducts(p,protons_all,protons_charm,protons_had,lambda,cBaryonDec);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the upsilons
      for (const Particle& ups : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==300553)) {
        Particles protons[3],lambda;
        findDecayProducts(ups,protons[0],protons[1],protons[2],lambda);
        LorentzTransform boost;
        if (ups.p3().mod() > 1*MeV) {
          boost = LorentzTransform::mkFrameTransformFromBeta(ups.mom().betaVec());
        }
        for (const Particle& p : protons[0]) {
          double pcm = boost.transform(p.mom()).p3().mod();
          _h_spect[0]->fill(pcm);
        }
        for (const Particle& p : lambda) {
          double pcm = boost.transform(p.mom()).p3().mod();
          _h_spect[1]->fill(pcm);
        }
        for (unsigned int ix=0; ix<2; ++ix) {
          for (const Particle& p : protons[ix+1]) {
            double pcm = boost.transform(p.mom()).p3().mod();
            _h_flav[ix]->fill(pcm);
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_spect, 1.0, false);
      normalize(_h_flav , 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_spect[2], _h_flav[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1992_I334962);

}
