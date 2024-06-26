// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> Xs l+ l-
  class BELLE_2005_I679052 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2005_I679052);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==521 || Cuts::abspid==511), "UFS");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix],1,1,1+ix);
      }
      book(_c,"TMP/Nb");
    }

    void findDecayProducts(bool& charm, const Particle& mother,
                           unsigned int& nK0, unsigned int& nKp,
			   unsigned int& nKm, Particles & lp, Particles & lm) {
      for (const Particle & p : mother.children()) {
        int id = p.pid();
        if (PID::isCharmHadron(p.pid())) charm = true;
        else if (id == PID::POSITRON || id == PID::ANTIMUON) lp.push_back(p);
        else if (id == PID::ELECTRON || id == PID::MUON    ) lm.push_back(p);
        else if (id == PID::KPLUS )    ++nKp;
        else if (id == PID::KMINUS )    ++nKm;
        else if (id == PID::K0S)        ++nK0;
        else if (id == PID::PI0 || id == PID::PIPLUS || id == PID::PIMINUS) {
          continue;
        }
        else if (!p.children().empty() ) {
          findDecayProducts(charm, p, nK0, nKp, nKm,lp,lm);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over bottoms
      for (const Particle& bottom : apply<UnstableParticles>(event, "UFS").particles()) {
       	// remove mixing entries etc
        if (bottom.children()[0].abspid()==bottom.abspid()) continue;
        _c->fill();
        bool charm = false;
        Particles lp,lm;
        unsigned int nK0(0),nKp(0),nKm(0);
        findDecayProducts(charm,bottom, nK0, nKp, nKm,lp,lm);
        if (charm) continue;
        unsigned int nk = nKp-nKm+nK0;
        if ( nk % 2 == 0) continue;
        if (lp.size()!=1 || lm.size()!=1 || lp[0].pid()!=-lm[0].pid()) continue;
        FourMomentum qq = lp[0].mom()+lm[0].mom();
        double q2 = qq.mass2();
        if (q2<0.04) continue;
        _h[0]->fill((bottom.mom()-qq).mass());
        _h[1]->fill(q2);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, 0.5e6/ *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c;
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2005_I679052);

}
