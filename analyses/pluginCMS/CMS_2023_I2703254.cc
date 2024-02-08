// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief ttbb in lepton+jets at 13 TeV
  class CMS_2023_I2703254 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2023_I2703254);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      Cut eta_cut = (Cuts::abseta < 5.0);
      const FinalState fs(eta_cut);

      PromptFinalState photons      (eta_cut && Cuts::abspid == PID::PHOTON, TauDecaysAs::PROMPT, MuDecaysAs::PROMPT);
      PromptFinalState electrons    (eta_cut && Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT, MuDecaysAs::PROMPT);
      PromptFinalState muons        (eta_cut && Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT, MuDecaysAs::PROMPT);

      Cut electron_cut     = (Cuts::abseta < 2.4) && (Cuts::pT > 29*GeV);
      Cut muon_cut         = (Cuts::abseta < 2.4) && (Cuts::pT > 26*GeV);
      Cut vetoelectron_cut = (Cuts::abseta < 2.5) && (Cuts::pT > 15*GeV);
      Cut vetomuon_cut     = (Cuts::abseta < 2.4) && (Cuts::pT > 15*GeV);

      LeptonFinder dressedelectrons    (electrons, photons, 0.1, electron_cut);
      LeptonFinder dressedmuons        (muons,     photons, 0.1, muon_cut);
      LeptonFinder dressedvetoelectrons(electrons, photons, 0.1, vetoelectron_cut);
      LeptonFinder dressedvetomuons    (muons,     photons, 0.1, vetomuon_cut);

      declare(dressedelectrons,     "elecs");
      declare(dressedmuons,         "muons");
      declare(dressedvetoelectrons, "vetoelecs");
      declare(dressedvetomuons,     "vetomuons");

      FastJets jetfs(fs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::NONE);
      declare(jetfs, "jets");

      // Book histograms
      book(_xsec, 1, 1, 1);
      book(_h["average_DR_6j_4b"],            11, 1, 1);
      book(_h["bjet3_abseta_5j_3b"],          12, 1, 1);
      book(_h["bjet3_abseta_6j_4b"],          13, 1, 1);
      book(_h["bjet3_pt_5j_3b"],              14, 1, 1);
      book(_h["bjet3_pt_6j_4b"],              15, 1, 1);
      book(_h["bjet4_abseta_6j_4b"],          16, 1, 1);
      book(_h["bjet4_pt_6j_4b"],              17, 1, 1);
      book(_h["bJet_Ht_5j_3b"],               18, 1, 1);
      book(_h["bJet_Ht_6j_4b"],               19, 1, 1);
      book(_h["dPhi_lj_bsoft_6j_3b_3lj"],     20, 1, 1);
      book(_h["dPhi_lj_bsoft_7j_4b_3lj"],     21, 1, 1);
      book(_h["extra_jet1_abseta_6j_4b"],     22, 1, 1);
      book(_h["extra_jet1_pt_6j_4b"],         23, 1, 1);
      book(_h["extra_jet2_abseta_6j_4b"],     24, 1, 1);
      book(_h["extra_jet2_pt_6j_4b"],         25, 1, 1);
      book(_h["extra_jet_abseta_6j_4b"],      26, 1, 1);
      book(_h["extra_jet_DR_6j_4b"],          27, 1, 1);
      book(_h["extra_jet_M_6j_4b"],           28, 1, 1);
      book(_h["extra_jet_pt_6j_4b"],          29, 1, 1);
      book(_h["extra_lightJet_pt_6j_3b_3lj"], 30, 1, 1);
      book(_h["extra_lightJet_pt_7j_4b_3lj"], 31, 1, 1);
      book(_h["jet_Ht_5j_3b"],                32, 1, 1);
      book(_h["jet_Ht_6j_4b"],                33, 1, 1);
      book(_h["largest_Mbb_6j_4b"],           34, 1, 1);
      book(_h["lightJet_Ht_6j_3b_3lj"],       35, 1, 1);
      book(_h["lightJet_Ht_7j_4b_3lj"],       36, 1, 1);
      book(_h["Njets_5j_3b"],                 37, 1, 1);
      book(_h["Njets_6j_4b"],                 38, 1, 1);
      book(_h["n_M_tagged_jets_5j_3b"],       39, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges.empty())  _edges = _xsec->xEdges();

      Particles leptons;
      for (auto &lep : apply<LeptonFinder>(event, "elecs").dressedLeptons()) { leptons.push_back(lep); }
      for (auto &lep : apply<LeptonFinder>(event, "muons").dressedLeptons()) { leptons.push_back(lep); }

      DressedLeptons Electrons     = apply<LeptonFinder>(event, "elecs").dressedLeptons();
      DressedLeptons Muons         = apply<LeptonFinder>(event, "muons").dressedLeptons();
      DressedLeptons vetoElectrons = apply<LeptonFinder>(event, "vetoelecs").dressedLeptons();
      DressedLeptons vetoMuons     = apply<LeptonFinder>(event, "vetomuons").dressedLeptons();


      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      const Jets nonisojets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.4);
      const Jets jets = discardIfAnyDeltaRLess(nonisojets, leptons, 0.4);

      Jets bjets;
      Jets lfjets;

      for (const Jet& jet : jets) {
        if (jet.bTagged())  bjets += jet;
        else lfjets += jet;
      }

      size_t njets   = jets.size();
      size_t nbjets  = bjets.size();
      size_t nlfjets = lfjets.size();

      bool pass_ljets   = (Electrons.size() == 1 && vetoElectrons.size() == 1 && vetoMuons.size()     == 0)
                           || (Muons.size() == 1 && vetoMuons.size()     == 1 && vetoElectrons.size() == 0);
      if (!(pass_ljets))            vetoEvent;
      if (nbjets < 3 || njets < 5)  vetoEvent;


      double ht_jets   = sum(  jets, Kin::pT, 0.0);
      double ht_bjets  = sum( bjets, Kin::pT, 0.0);
      double ht_lfjets = sum(lfjets, Kin::pT, 0.0);

      size_t ind1{0}, ind2{0}; double mindr = 999.;
      double maxmbb  =   0.;
      double sum_dr = 0.0;
      size_t sum_n_dr = 0;
      for (size_t i = 0; i < bjets.size(); ++i) {
        for (size_t j = i + 1; j < bjets.size(); ++j) {
          double dr = deltaR(bjets[i], bjets[j]);
          double mbb = (bjets[i].momentum() + bjets[j].momentum()).mass();

          if (dr < mindr) {
            ind1 = i;
            ind2 = j;
            mindr = dr;
          }
          sum_dr += dr;
          sum_n_dr += 1;

          if (mbb > maxmbb) maxmbb = mbb;
        }
      }

      FourMomentum bb_closest = bjets[ind1].momentum() + bjets[ind2].momentum();
      double dr_closest = deltaR(bjets[ind1], bjets[ind2]);

      // 6j3b3l and 7j4b3l:
      double gen_w_mass = 79.6;
      size_t lfind1{0}, lfind2{0}; double minwmasdiff = 9999.;
      size_t lfind = 0;
      double phi_first_ljet_extra_bjet_soft = 0.;
      if (nlfjets >= 3) {
        for (size_t i = 0; i < nlfjets; ++i) {
          for (size_t j = i + 1; j < nlfjets; ++j) {
            double wmasdiff = abs(gen_w_mass - (lfjets[i].momentum()+lfjets[j].momentum()).mass());
            if (wmasdiff < minwmasdiff) {
              lfind1 = i;
              lfind2 = j;
              minwmasdiff = wmasdiff;
            }
          }
        }
        for (size_t i = 0; i < 3; ++i){
            if (!(( i == lfind1) || ( i == lfind2))){
            lfind = i;
            break;
            }
        }
        phi_first_ljet_extra_bjet_soft = deltaPhi(lfjets[lfind].momentum(),bjets[nbjets - 1].momentum());
      }

      // Fill histograms
      // 5j3b category
      if (njets >= 5 && nbjets >= 3) {
        _xsec->fill(_edges[3]);
        _h["Njets_5j_3b"]                   ->fill(njets);
        _h["n_M_tagged_jets_5j_3b"]         ->fill(nbjets);
        _h["bjet3_pt_5j_3b"]                ->fill(bjets[2].pT()/GeV);
        _h["bjet3_abseta_5j_3b"]            ->fill(bjets[2].abseta()/GeV);
        _h["jet_Ht_5j_3b"]                  ->fill(ht_jets/GeV);
        _h["bJet_Ht_5j_3b"]                 ->fill(ht_bjets/GeV);
      }

      // 6j4b category
      if (njets >= 6 && nbjets >= 4) {
        _xsec->fill(_edges[1]);
        _h["Njets_6j_4b"]                   ->fill(njets);
        _h["bjet3_pt_6j_4b"]                ->fill(bjets[2].pT()/GeV);
        _h["bjet3_abseta_6j_4b"]            ->fill(bjets[2].abseta()/GeV);
        _h["bjet4_pt_6j_4b"]                ->fill(bjets[3].pT()/GeV);
        _h["bjet4_abseta_6j_4b"]            ->fill(bjets[3].abseta()/GeV);
        _h["jet_Ht_6j_4b"]                  ->fill(ht_jets/GeV);
        _h["bJet_Ht_6j_4b"]                 ->fill(ht_bjets/GeV);
        _h["average_DR_6j_4b"]              ->fill(sum_dr/sum_n_dr);
        _h["largest_Mbb_6j_4b"]             ->fill(maxmbb);
        _h["extra_jet_DR_6j_4b"]            ->fill(dr_closest);
        _h["extra_jet_M_6j_4b"]             ->fill(bb_closest.mass()/GeV);
        _h["extra_jet_pt_6j_4b"]            ->fill(bb_closest.pT()/GeV);
        _h["extra_jet_abseta_6j_4b"]        ->fill(bb_closest.abseta()/GeV);
        _h["extra_jet1_pt_6j_4b"]           ->fill(bjets[ind1].pT()/GeV);
        _h["extra_jet1_abseta_6j_4b"]       ->fill(bjets[ind1].abseta()/GeV);
        _h["extra_jet2_pt_6j_4b"]           ->fill(bjets[ind2].pT()/GeV);
        _h["extra_jet2_abseta_6j_4b"]       ->fill(bjets[ind2].abseta()/GeV);
      }
      // 6j3b3l category
      if (njets >= 6 && nbjets >= 3 && nlfjets >= 3) {
        _xsec->fill(_edges[2]);
        _h["extra_lightJet_pt_6j_3b_3lj"]   ->fill(lfjets[lfind].pT()/GeV);
        _h["dPhi_lj_bsoft_6j_3b_3lj"]       ->fill(phi_first_ljet_extra_bjet_soft);
        _h["lightJet_Ht_6j_3b_3lj"]         ->fill(ht_lfjets/GeV);
      }

      // 7j4b3l category
      if (njets >= 7 && nbjets >= 4 && nlfjets >= 3) {
        _xsec->fill(_edges[0]);
        _h["extra_lightJet_pt_7j_4b_3lj"]   ->fill(lfjets[lfind].pT()/GeV);
        _h["dPhi_lj_bsoft_7j_4b_3lj"]       ->fill(phi_first_ljet_extra_bjet_soft);
        _h["lightJet_Ht_7j_4b_3lj"]         ->fill(ht_lfjets/GeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection() / femtobarn / sumOfWeights();

      // move the overflow into the last "true" bin of the distribution _just_ for the following
      const std::list<string> overflows =
      {
        "Njets_5j_3b",
        "bJet_Ht_5j_3b",
        "bjet3_pt_5j_3b",
        "jet_Ht_5j_3b",
        "n_M_tagged_jets_5j_3b",
        "extra_lightJet_pt_6j_3b_3lj",
        "lightJet_Ht_6j_3b_3lj",
        "Njets_6j_4b",
        "bJet_Ht_6j_4b",
        "bjet3_pt_6j_4b",
        "bjet4_pt_6j_4b",
        "extra_jet1_pt_6j_4b",
        "extra_jet2_pt_6j_4b",
        "extra_jet_M_6j_4b",
        "extra_jet_abseta_6j_4b",
        "extra_jet_pt_6j_4b",
        "jet_Ht_6j_4b",
        "largest_Mbb_6j_4b",
        "extra_lightJet_pt_7j_4b_3lj",
        "lightJet_Ht_7j_4b_3lj"
      };
      for (const string& of : overflows) {
        auto& of_bin = _h[of]->bin(_h[of]->numBins()+1);
        _h[of]->bin(_h[of]->numBins()) += of_bin;
        of_bin.reset();
      }

      scale(_xsec, sf);
      normalize(_h);

    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    BinnedHistoPtr<string> _xsec;
    vector<string> _edges;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(CMS_2023_I2703254);

}
