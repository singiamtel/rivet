// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// @brief ATLAS 13 TeV Z(->ll)yy analysis
  class ATLAS_2022_I2593322 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2022_I2593322);

    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Prompt photons
      const Cut photoncut = Cuts::abspid == PID::PHOTON && Cuts::pT > 20*GeV && Cuts::abseta < 2.37;
      PromptFinalState photon_fs(photoncut);
      declare(photon_fs, "Photons");

      // Prompt leptons
      const PromptFinalState bareelectron_fs = Cuts::abspid == PID::ELECTRON;
      const PromptFinalState baremuon_fs = Cuts::abspid == PID::MUON;

      // Dressed leptons
      const IdentifiedFinalState allphoton_fs(PID::PHOTON); // photons used for lepton dressing
      const Cut leptoncut = Cuts::pT > 20*GeV && Cuts::abseta < 2.47;
      const DressedLeptons dressedelectron_fs(allphoton_fs, bareelectron_fs, 0.1, leptoncut, true); // use *all* photons for lepton dressing
      const DressedLeptons dressedmuon_fs(allphoton_fs, baremuon_fs, 0.1, leptoncut, true); // use *all* photons for lepton dressing

      declare(dressedelectron_fs, "Electrons");
      declare(dressedmuon_fs, "Muons");

      IdentifiedFinalState neutrinos;
      neutrinos.acceptNeutrinos();

      // FS for photon isolation
      FinalState all_fs;
      VetoedFinalState veto_fs(all_fs);
      veto_fs.addVetoOnThisFinalState(photon_fs);
      veto_fs.addVetoOnThisFinalState(dressedmuon_fs);
      veto_fs.addVetoOnThisFinalState(neutrinos);
      declare(veto_fs, "vetoFS");

      // book histograms
      book(_h["y1_pt"], 2, 1, 1);
      book(_h["y2_pt"], 3, 1, 1);
      book(_h["ll_pt"], 4, 1, 1);
      book(_h["llyy_pt"], 5, 1, 1);
      book(_h["yy_m"], 6, 1, 1);
      book(_h["llyy_m"], 7, 1, 1);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      vector<DressedLepton> electrons = apply<DressedLeptons>(event, "Electrons").dressedLeptons();
      vector<DressedLepton> muons = apply<DressedLeptons>(event, "Muons").dressedLeptons();
      const Particles& photons = apply<PromptFinalState>(event, "Photons").particlesByPt();

      if ( (electrons.size() < 2 && muons.size() < 2) ) vetoEvent;
      if ( photons.size() < 2 ) vetoEvent;

      DressedLepton *lep_1, *lep_2;

      if (muons.size() >= 2){
        if (muons[0].pT()/GeV < 30) vetoEvent;
        lep_1 = &muons[0];
        lep_2 = &muons[1];
      }
      else{
        if (electrons[0].pT()/GeV < 30) vetoEvent;
        lep_1 = &electrons[0];
        lep_2 = &electrons[1];
      }

      if ( (lep_1->charge() == lep_2->charge()) || (lep_1->abspid() != lep_2->abspid()) ) vetoEvent;

      Particles veto_particles = apply<VetoedFinalState>(event, "vetoFS").particles();
      Particles selPhotons;
      for (size_t i = 0; i < photons.size(); ++i){
        if ( deltaR(photons[i], *lep_1) < 0.4 )  continue;
        if ( deltaR(photons[i], *lep_2) < 0.4 )  continue;
        double coneEnergy = 0;
        for (const Particle &p : veto_particles){
          if ( deltaR(photons[i], p) < 0.2 )  coneEnergy += p.Et();
        }
        if (coneEnergy/photons[i].Et() > 0.07) continue;
        selPhotons.push_back(photons[i]);
      }

      if (selPhotons.size() < 2) vetoEvent;

      if (deltaR(selPhotons[0], selPhotons[1]) < 0.4) vetoEvent;

      if ( (lep_1->momentum() + lep_2->momentum()).mass()/GeV < 40 ) vetoEvent;

      FourMomentum ll_p4 = lep_1->mom()+lep_2->mom();
      const double m_ll = ll_p4.mass();
      const bool lly0 = (ll_p4+selPhotons[0].mom()).mass() > (ll_p4+selPhotons[1].mom()).mass();
      const double m_ll_y = lly0? (ll_p4+selPhotons[1].mom()).mass() : (ll_p4+selPhotons[0].mom()).mass();

      if ( (m_ll + m_ll_y) < 182*GeV)  vetoEvent;

      FourMomentum yy_p4 = selPhotons[0].momentum()+selPhotons[1].momentum();

      _h["y1_pt"]->fill(selPhotons[0].pt()/GeV);
      _h["y2_pt"]->fill(selPhotons[1].pt()/GeV);
      _h["ll_pt"]->fill(ll_p4.pt()/GeV);
      _h["llyy_pt"]->fill( (ll_p4+yy_p4).pt()/GeV );
      _h["yy_m"]->fill(yy_p4.mass()/GeV);
      _h["llyy_m"]->fill( (ll_p4+yy_p4).mass()/GeV );

    }

    void finalize() {

      const double sf = crossSection()/femtobarn/sumW();
      scale(_h, sf);

    }

    ///@}

    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2022_I2593322);

}
