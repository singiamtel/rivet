// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// D0 direct photon pair production
  class D0_2010_I846997 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2010_I846997);


    void init() {
      FinalState fs;
      declare(fs, "FS");

      IdentifiedFinalState ifs(Cuts::abseta < 0.9 && Cuts::pT > 20*GeV);
      ifs.acceptId(PID::PHOTON);
      declare(ifs, "IFS");

      book(_h_M ,1, 1, 1);
      book(_h_pT ,2, 1, 1);
      book(_h_dPhi ,3, 1, 1);
      book(_h_costheta ,4, 1, 1);

      book(_h_pT_M, {30., 50., 80., 350.}, {"d05-x01-y01", "d08-x01-y01", "d11-x01-y01"});
      book(_h_dPhi_M, {30., 50., 80., 350.}, {"d06-x01-y01", "d09-x01-y01", "d12-x01-y01"});
      book(_h_costheta_M, {30., 50., 80., 350.}, {"d07-x01-y01", "d10-x01-y01", "d13-x01-y01"});
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Particles photons = apply<IdentifiedFinalState>(event, "IFS").particlesByPt();
      if (photons.size() < 2 ||
          (photons[0].pT() < 21.0*GeV)) {
        vetoEvent;
      }

      // Isolate photons with ET_sum in cone
      Particles isolated_photons;
      Particles fs = apply<FinalState>(event, "FS").particles();
      for (const Particle& photon : photons) {
        double eta_P = photon.eta();
        double phi_P = photon.phi();
        double Etsum=0.0;
        for (const Particle& p : fs) {
          if (HepMCUtils::uniqueId(p.genParticle()) != HepMCUtils::uniqueId(photon.genParticle()) &&
              deltaR(eta_P, phi_P, p.eta(), p.phi()) < 0.4) {
            Etsum += p.Et();
          }
        }
        if (Etsum < 2.5*GeV) {
          isolated_photons.push_back(photon);
        }
      }

      if (isolated_photons.size() != 2) {
        vetoEvent;
      }
      std::sort(isolated_photons.begin(), isolated_photons.end(), cmpMomByPt);

      FourMomentum y1=isolated_photons[0].momentum();
      FourMomentum y2=isolated_photons[1].momentum();
      if (deltaR(y1, y2)<0.4)  vetoEvent;

      FourMomentum yy=y1+y2;
      double Myy = yy.mass()/GeV;
      if (Myy<30.0 || Myy>350.0)  vetoEvent;

      double pTyy = yy.pT()/GeV;
      if (Myy<pTyy)  vetoEvent;

      double dPhiyy = mapAngle0ToPi(y1.phi()-y2.phi());
      if (dPhiyy<0.5*M_PI)  vetoEvent;

      double costhetayy = fabs(tanh((y1.eta()-y2.eta())/2.0));

      _h_M->fill(Myy);
      _h_pT->fill(pTyy);
      _h_dPhi->fill(dPhiyy);
      _h_costheta->fill(costhetayy);

      _h_pT_M->fill(Myy, pTyy);
      _h_dPhi_M->fill(Myy, dPhiyy);
      _h_costheta_M->fill(Myy, costhetayy);
    }


    void finalize() {

      scale(_h_M, crossSection()/picobarn/sumOfWeights());
      scale(_h_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_dPhi, crossSection()/picobarn/sumOfWeights());
      scale(_h_costheta, crossSection()/picobarn/sumOfWeights());

      scale(_h_pT_M, crossSection()/picobarn/sumOfWeights());
      scale(_h_dPhi_M, crossSection()/picobarn/sumOfWeights());
      scale(_h_costheta_M, crossSection()/picobarn/sumOfWeights());

      divByGroupWidth({_h_pT_M, _h_dPhi_M, _h_costheta_M});
    }


  private:

    Histo1DPtr _h_M;
    Histo1DPtr _h_pT;
    Histo1DPtr _h_dPhi;
    Histo1DPtr _h_costheta;
    Histo1DGroupPtr _h_pT_M, _h_dPhi_M, _h_costheta_M;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2010_I846997, D0_2010_S8570965);

}
