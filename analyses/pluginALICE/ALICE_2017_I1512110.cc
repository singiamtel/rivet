//-*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  class ALICE_2017_I1512110 : public Analysis {
  public:

    /// Constructor
    ALICE_2017_I1512110()
      : Analysis("ALICE_2017_I1512110"),
        _rapmax(0.8)
    {    }


    void init() {

      const UnstableParticles ufs(Cuts::absrap < _rapmax);
      declare(ufs, "UFS");

      book(_h_pi0, 3,1,1);
      book(_h_eta, 4,1,1);
      book(_h_etaToPion, 5,1,1);

      // temporary plots with the binning of _h_etaToPion
      // to construct the eta/pi0 ratio in the end
      book(_temp_h_pion, "TMP/h_pion",refData(5,1,1));
      book(_temp_h_eta , "TMP/h_eta", refData(5,1,1));
    }


    void analyze(const Event& event) {

      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles()) {

        if (p.pid() == 111) {
          // neutral pion; ALICE corrects for pi0 feed-down
          if ( !(p.hasAncestorWith(Cuts::pid == 310)  || p.hasAncestorWith(Cuts::pid == 130)   || // K0_s, K0_l
                 p.hasAncestorWith(Cuts::pid == 321)  || p.hasAncestorWith(Cuts::pid == -321)  || // K+,K-
                 p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) || // Lambda, Anti-Lambda
                 p.hasAncestorWith(Cuts::pid == 3212) || p.hasAncestorWith(Cuts::pid == -3212) || // Sigma0
                 p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) || // Sigmas
                 p.hasAncestorWith(Cuts::pid == 3112) || p.hasAncestorWith(Cuts::pid == -3112) || // Sigmas
                 p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) || // Cascades
                 p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) )) // Cascades
            {
              _h_pi0->fill(p.pT()/GeV, 1.0 /(TWOPI*p.pT()/GeV*2*_rapmax));
              _temp_h_pion->fill(p.pT()/GeV);
            }
        }
        else if (p.pid() == 221){
          // eta meson
          _h_eta->fill(p.pT()/GeV, 1.0 /(TWOPI*p.pT()/GeV*2*_rapmax));
          _temp_h_eta->fill(p.pT()/GeV);
        }
      }
    }


    void finalize() {

      scale(_h_pi0, crossSection()/picobarn/sumOfWeights());
      scale(_h_eta, crossSection()/picobarn/sumOfWeights());
      divide(_temp_h_eta, _temp_h_pion, _h_etaToPion);

    }


  private:

    double _rapmax;
    Histo1DPtr _h_pi0;
    Histo1DPtr _h_eta;
    Histo1DPtr _temp_h_pion;
    Histo1DPtr _temp_h_eta;
    Estimate1DPtr _h_etaToPion;

  };


  RIVET_DECLARE_PLUGIN(ALICE_2017_I1512110);

}
