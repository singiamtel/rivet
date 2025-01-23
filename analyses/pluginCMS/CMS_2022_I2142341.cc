// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/RivetYODA.hh"

namespace Rivet {

  /// @brief H->yy analysis at 13 TeV
  class CMS_2022_I2142341 : public Analysis {
   public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2022_I2142341);

    void init() {
      //---All particle final state
      FinalState fs;
      declare(fs, "FS");

      //---Photons
      FinalState fs_photons(Cuts::abspid == PID::PHOTON);
      declare(fs_photons, "FS_PHOTONS");

      //---Jets
      FastJets fs_jets(fs, JetAlg::ANTIKT, 0.4);
      declare(fs_jets, "JETS");

      //---Histograms
      book(_h["pt_gg"], "d01-x01-y01");
      book(_h["njets_eta2p5"], "d03-x01-y01");
      book(_h["cos_theta_star"], "d05-x01-y01");
      book(_h["rapidity_gg"], "d07-x01-y01");
      book(_h["jet_pt_eta2p5"], "d09-x01-y01");
      book(_h["jet_rapidity_eta2p5"], "d11-x01-y01");
      book(_h["deltaphijj"], "d23-x01-y01");
      book(_h["deltaetajj"], "d29-x01-y01");
    }

    // cos theta star angle in the Collins Soper frame
    double getCosThetaStar_CS(const FourMomentum& h1, const FourMomentum& h2) {
      FourMomentum hh = h1 + h2;
      LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(hh.betaVec());
      FourMomentum h1_boosted = boost.transform(h1);
      return abs(cos(h1_boosted.theta()));
    }

    void analyze(const Event& event) {
      Particles photons = apply<FinalState>(event, "FS_PHOTONS").particlesByPt();

      if (photons.size() < 2) vetoEvent;

      //---Isolate photons with ET_sum in cone
      Particles isolated_photons;
      Particles fs = apply<FinalState>(event, "FS").particlesByPt();
      for (const Particle& photon : photons) {
        FourMomentum mom_in_cone;
        for (const Particle& particle : fs) {
          if (deltaR(photon, particle) < 0.3 && photon.p() != particle.p())
            mom_in_cone += particle.momentum();
        }
        if (mom_in_cone.pt() < 10 * GeV && photon.abseta() < 2.5)
          isolated_photons.push_back(photon);
      }

      if (isolated_photons.size() < 2)
        vetoEvent;

      //---kinematic photon selection
      FourMomentum mom_PP = isolated_photons[0].mom() + isolated_photons[1].mom();
      if (isolated_photons[0].pT() < mom_PP.mass() / 3. || isolated_photons[1].pT() < mom_PP.mass() / 4.)
        vetoEvent;

      _h["pt_gg"]->fill(mom_PP.pt() / GeV);
      _h["rapidity_gg"]->fill(abs(mom_PP.rapidity()));
      _h["cos_theta_star"]->fill(getCosThetaStar_CS(isolated_photons[0].mom(), isolated_photons[1].mom()));

      //---jets
      auto jets_eta4p7 = apply<FastJets>(event, "JETS").jetsByPt(Cuts::abseta < 4.7 && Cuts::pt > 30 * GeV);
      Jets jets_eta4p7_npart;
      for (const Jet& j : jets_eta4p7) {
        if (j.constituents().size() > 5) jets_eta4p7_npart.push_back(j);
      }
      Particles selected_photons;
      selected_photons.push_back(isolated_photons[0]);
      selected_photons.push_back(isolated_photons[1]);
      idiscardIfAnyDeltaRLess(jets_eta4p7_npart, selected_photons, 0.4);
      if (jets_eta4p7_npart.size() > 1) {
        _h["deltaphijj"]->fill(deltaPhi(jets_eta4p7_npart[0], jets_eta4p7_npart[1]));
        _h["deltaetajj"]->fill(deltaEta(jets_eta4p7_npart[0], jets_eta4p7_npart[1]));
      }

      Jets jets_eta2p5;
      for (const Jet& j : jets_eta4p7_npart) {
        if (j.abseta() < 2.5) jets_eta2p5.push_back(j);
      }
      _h["njets_eta2p5"]->fill(jets_eta2p5.size());
      if (jets_eta2p5.size() > 0) {
        _h["jet_pt_eta2p5"]->fill(jets_eta2p5[0].pt() / GeV);
        _h["jet_rapidity_eta2p5"]->fill(abs(jets_eta2p5[0].rapidity()));
      }
    }

    void finalize() {
      scale(_h, crossSection() / femtobarn * BR / sumOfWeights());
    }

   private:

    map<string, Histo1DPtr> _h;
    const double BR = 0.00227;
  };

  RIVET_DECLARE_PLUGIN(CMS_2022_I2142341);

}  // namespace Rivet
