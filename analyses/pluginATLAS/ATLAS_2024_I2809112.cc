// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"

namespace Rivet {


class ATLAS_2024_I2809112 : public Analysis {
 public:

  /// @brief ttbb production at 13 TeV
  RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2024_I2809112);

  void init() {

    // Photons
    PromptFinalState photons(Cuts::abspid == PID::PHOTON);

    // Kinematic cuts for leptons
    const Cut cuts_lep = Cuts::pT > 25*GeV && Cuts::abseta < 2.5;

    // Muons
    PromptFinalState bare_mu(Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT);
    LeptonFinder all_dressed_mu(bare_mu, photons, 0.1, cuts_lep);
    declare(all_dressed_mu, "DressedMuons");

    // Electrons
    PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);
    LeptonFinder all_dressed_el(bare_el, photons, 0.1, cuts_lep);
    declare(all_dressed_el, "DressedElectrons");

    //Jet forming
    const InvisibleFinalState neutrinos(OnlyPrompt::YES, TauDecaysAs::PROMPT);

    VetoedFinalState vfs(FinalState(Cuts::abseta < 5.0));
    vfs.addVetoOnThisFinalState(all_dressed_el);
    vfs.addVetoOnThisFinalState(all_dressed_mu);
    vfs.addVetoOnThisFinalState(neutrinos);

    FastJets jet(vfs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::ALL);
    declare(jet, "Jets");

    const FinalState fs(Cuts::abseta < 5.);
    declare(MissingMomentum(fs), "ETmiss");

    // Njets
    book(_d["nbjets_emu2b"],     71, 1, 1);
    book(_d["nlightjets_emu3b"], 95, 1, 1);
    book(_d["nlightjets_emu4b"], 134, 1, 1);

    // HT
    book(_h["ht_had_emu3b"],  72, 1, 1);
    book(_h["ht_all_emu3b"],  73, 1, 1);
    book(_h["ht_had_emu4b"], 102, 1, 1);
    book(_h["ht_all_emu4b"], 103, 1, 1);

    // b-jet pTs
    book(_h["lead_bjet_pt_emu3b"],     76, 1, 1);
    book(_h["sublead_bjet_pt_emu3b"],  78, 1, 1);
    book(_h["third_bjet_pt_emu3b"],    80, 1, 1);
    book(_h["lead_bjet_pt_emu4b"],    106, 1, 1);
    book(_h["sublead_bjet_pt_emu4b"], 108, 1, 1);
    book(_h["third_bjet_pt_emu4b"],   110, 1, 1);
    book(_h["fourth_bjet_pt_emu4b"],  112, 1, 1);
    // bjet eta
    book(_h["lead_bjet_eta_emu3b"],     82, 1, 1);
    book(_h["sublead_bjet_eta_emu3b"],  84, 1, 1);
    book(_h["third_bjet_eta_emu3b"],    86, 1, 1);
    book(_h["lead_bjet_eta_emu4b"],    114, 1, 1);
    book(_h["sublead_bjet_eta_emu4b"], 116, 1, 1);
    book(_h["third_bjet_eta_emu4b"],   118, 1, 1);
    book(_h["fourth_bjet_eta_emu4b"],  120, 1, 1);

    //bjets from top with limited accuracy
    book(_h["lead_bjetfromtop_pt_emu3b"],      77, 1, 1);
    book(_h["sublead_bjetfromtop_pt_emu3b"],   79, 1, 1);
    book(_h["lead_bjetfromtop_pt_emu4b"],     107, 1, 1);
    book(_h["sublead_bjetfromtop_pt_emu4b"],  109, 1, 1);
    book(_h["lead_bjetfromtop_eta_emu3b"],     83, 1, 1);
    book(_h["sublead_bjetfromtop_eta_emu3b"],  85, 1, 1);
    book(_h["lead_bjetfromtop_eta_emu4b"],    115, 1, 1);
    book(_h["sublead_bjetfromtop_eta_emu4b"], 117, 1, 1);

    //additional bjets
    book(_h["lead_bjetfromQCD_pt_emu3b"],      81, 1, 1);
    book(_h["lead_bjetfromQCD_eta_emu3b"],     87, 1, 1);
    book(_h["lead_bjetfromQCD_pt_emu4b"],     111, 1, 1);
    book(_h["sublead_bjetfromQCD_pt_emu4b"],  113, 1, 1);
    book(_h["lead_bjetfromQCD_eta_emu4b"],    119, 1, 1);
    book(_h["sublead_bjetfromQCD_eta_emu4b"], 121, 1, 1);

    // leading bb pair
    book(_h["m_bb_leading_emu3b"],   88, 1, 1);
    book(_h["pt_bb_leading_emu3b"],  89, 1, 1);
    book(_h["dR_bb_leading_emu3b"],  94, 1, 1);
    book(_h["m_bb_leading_emu4b"],  122, 1, 1);
    book(_h["pt_bb_leading_emu4b"], 123, 1, 1);
    book(_h["dR_bb_leading_emu4b"], 133, 1, 1);

    // closest bb pair
    book(_h["m_bb_closest_emu4b"],  128, 1, 1);
    book(_h["pt_bb_closest_emu4b"], 129, 1, 1);
    book(_h["dR_bb_closest_emu4b"], 132, 1, 1);

    // bb pair from top
    book(_h["m_bbfromtop_emu3b"],   90, 1, 1);
    book(_h["pt_bbfromtop_emu3b"],  91, 1, 1);
    book(_h["m_bbfromtop_emu4b"],  124, 1, 1);
    book(_h["pt_bbfromtop_emu4b"], 125, 1, 1);

    book(_h["m_bbadd_emu4b"],  130, 1, 1);
    book(_h["pt_bbadd_emu4b"], 131, 1, 1);

    // light-jets
    book(_h["lead_lightjet_eta_emu3b"], 100, 1, 1);
    book(_h["lead_lightjet_pt_emu3b"],  101, 1, 1);
    book(_h["lead_lightjet_eta_emu4b"], 139, 1, 1);
    book(_h["lead_lightjet_pt_emu4b"],  140, 1, 1);

    // global variables
    book(_h["dR_avg_emu3b"],       74, 1, 1);
    book(_h["deta_jj_max_emu3b"],  75, 1, 1);
    book(_h["dR_avg_emu4b"],      104, 1, 1);
    book(_h["deta_jj_max_emu4b"], 105, 1, 1);

    book(_h["m_emubblead_emu3b"],    92, 1, 1);
    book(_h["m_emubbfromtop_emu3b"], 93, 1, 1);

    book(_h["dR_emubblead_thirdbjet_emu3b"],          96, 1, 1);
    book(_h["dR_emubbfromtop_leadaddb_emu3b"],        97, 1, 1);
    book(_h["dR_emubbfromtop_leadlightjet_emu3b"],    98, 1, 1);
    book(_h["ptdiff_leadaddbjet_leadlightjet_emu3b"], 99, 1, 1);

    book(_h["m_emubblead_emu4b"],    126, 1, 1);
    book(_h["m_emubbfromtop_emu4b"], 127, 1, 1);

    book(_h["dR_emubblead_thirdbjet_emu4b"],          135, 1, 1);
    book(_h["dR_emubbfromtop_leadaddb_emu4b"],        136, 1, 1);
    book(_h["dR_emubbfromtop_leadlightjet_emu4b"],    137, 1, 1);
    book(_h["ptdiff_leadaddbjet_leadlightjet_emu4b"], 138, 1, 1);

  }


  void analyze(const Event& event) {

    DressedLeptons muons = apply<LeptonFinder>(event, "DressedMuons").dressedLeptons();
    DressedLeptons electrons = apply<LeptonFinder>(event, "DressedElectrons").dressedLeptons();
    const double met = apply<MissingMomentum>(event, "ETmiss").missingPt();

    // lepton-jet overlap removal (note: muons are not included in jet finding)
    Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

    // Remove leptons within dR < 0.4 of a jet
    idiscardIfAnyDeltaRLess(muons, jets, 0.4);
    idiscardIfAnyDeltaRLess(electrons, jets, 0.4);

    Jets bjets, lightjets;
    for (const Jet& jet: jets) {
      if (jet.bTagged(Cuts::pT > 5.0*GeV))  bjets += jet;
      else                                  lightjets += jet;
    }

    // Evaluate basic event selection
    size_t nElecs = electrons.size(), nMuons = muons.size();
    const bool pass_emu = (nElecs == 1) && (nMuons ==1) &&
                          (electrons[0].pT() >= 28*GeV) && (muons[0].pT() >= 28*GeV) &&
                          (electrons[0].charge() * muons[0].charge() < 0.);


    // If we don't have exactly 2 leptons then veto the event
    size_t nbjets = bjets.size();
    if (!pass_emu || nbjets < 2) vetoEvent;

    _d["nbjets_emu2b"]->fill(nbjets>=5? "≥5"s : toString(nbjets));
    if (nbjets < 3) vetoEvent;

    DressedLeptons leptons = electrons + muons;
    const double hthad = sum(jets, Kin::pT, 0.0);
    const double htvis = sum(leptons, Kin::pT, hthad);
    const double htall = htvis + met;

    const FourMomentum jsum = bjets[0].mom() + bjets[1].mom();
    const double dr_leading = deltaR(bjets[0], bjets[1]);

    size_t idx1 = 0, idx2 = 1, npairs = 0;
    double dr_closest = 1e3, dr_avg = 0., maxdR_bb = 0.;
    for (size_t i = 0; i < bjets.size(); ++i) {
      for (size_t j = i+1; j < bjets.size(); ++j) {
        const double dr = deltaR(bjets[i], bjets[j]);
        dr_avg += dr;
        ++npairs;
        if (dr >= maxdR_bb)  maxdR_bb = dr;
        if (dr < dr_closest) {
          idx1 = i;
          idx2 = j;
          dr_closest = dr;
        }
      }
    }
    dr_avg /= npairs;

    const FourMomentum bb_closest = bjets[idx1].mom() + bjets[idx2].mom();
    const double mindR_l1b = drlb_min(leptons[0], bjets);
    const double mindR_l2b = drlb_min(leptons[1], bjets);

    double deta_max = 0.;
    for (size_t i = 0; i < jets.size(); ++i) {
      for (size_t j = i+1; j < jets.size(); ++j) {
        const double abs_eta = fabs(deltaEta(jets[i], jets[j]));
        if (abs_eta >= deta_max)  deta_max = abs_eta;
      }
    }

    const FourMomentum emubb = leptons[0].mom() + leptons[1].mom() + bjets[0].mom() + bjets[1].mom();
    const double mass_emubb = emubb.mass();
    const double dR_emubblead_thirdbjet = deltaR(emubb, bjets[2].mom());

    // permuted set of jets
    vector<Jets> bjet_sets;
    do {
      bjet_sets.push_back(bjets);
    } while (next_permutation(bjets.begin(), bjets.end()));

    //find jet pairs from top
    double m_bestweight_max = -9.;
    Jets m_bjets_bestweight;
    for (const Jets& bjet_vec : bjet_sets) {
    	if (bjet_vec.size() < 2) continue;
    	// calculate bjet from top weights using algorithm
      const double bfromtop_weight = bjetOriginWeight(leptons, bjet_vec, dr_closest,
                                                      maxdR_bb, mindR_l1b, mindR_l2b);
      if (bfromtop_weight > m_bestweight_max) {
        m_bestweight_max = bfromtop_weight;
        m_bjets_bestweight = bjet_vec;
      }
    }

    Jets bfromtop, bfromQCD;
    bfromtop += m_bjets_bestweight[0];
    bfromtop += m_bjets_bestweight[1];
    bfromQCD += m_bjets_bestweight[2];
    if (m_bjets_bestweight.size()>=4) bfromQCD += m_bjets_bestweight[3];

    isortByPt(bfromtop);
    isortByPt(bfromQCD);

    const FourMomentum jfromtopsum = bfromtop[0].mom() + bfromtop[1].mom();
    const FourMomentum dilepton = leptons[0].mom() + leptons[1].mom();
    const FourMomentum lljjfromtop = dilepton + jfromtopsum;
    const double mass_emubbfromtop = lljjfromtop.mass();
    const double dR_emubbfromtop_leadadd = deltaR(lljjfromtop, bfromQCD[0].mom());


    //nlight jets
    size_t nlightjets = lightjets.size();
    _d["nlightjets_emu3b"]->fill(nlightjets>=3? "≥3"s : toString(nlightjets));
    // b-jet pTs
    _h["lead_bjet_pt_emu3b"]->fill(bjets[0].pT()/GeV);
    _h["sublead_bjet_pt_emu3b"]->fill(bjets[1].pT()/GeV);
    _h["third_bjet_pt_emu3b"]->fill(bjets[2].pT()/GeV);
    // b-jet eta
    _h["lead_bjet_eta_emu3b"]->fill(fabs(bjets[0].eta()));
    _h["sublead_bjet_eta_emu3b"]->fill(bjets[1].abseta());
    _h["third_bjet_eta_emu3b"]->fill(bjets[2].abseta());

    // bjets assigned to top
    _h["lead_bjetfromtop_pt_emu3b"]->fill(bfromtop[0].pT()/GeV);
    _h["sublead_bjetfromtop_pt_emu3b"]->fill(bfromtop[1].pT()/GeV);
    _h["lead_bjetfromtop_eta_emu3b"]->fill(bfromtop[0].abseta());
    _h["sublead_bjetfromtop_eta_emu3b"]->fill(bfromtop[1].abseta());

    //additional bjets
    _h["lead_bjetfromQCD_pt_emu3b"]->fill(bfromQCD[0].pT()/GeV);
    _h["lead_bjetfromQCD_eta_emu3b"]->fill(bfromQCD[0].abseta());

    // HT
    _h["ht_had_emu3b"]->fill(hthad/GeV);
    _h["ht_all_emu3b"]->fill(htall/GeV);

    // leading bb pair
    _h["m_bb_leading_emu3b"]->fill(jsum.mass()/GeV);
    _h["pt_bb_leading_emu3b"]->fill(jsum.pT()/GeV);
    _h["dR_bb_leading_emu3b"]->fill(dr_leading);

    // b-jets pair assigned to top
    _h["m_bbfromtop_emu3b"]->fill(jfromtopsum.mass()/GeV);
    _h["pt_bbfromtop_emu3b"]->fill(jfromtopsum.pT()/GeV);

    // global
    _h["dR_avg_emu3b"]->fill(dr_avg);
    _h["deta_jj_max_emu3b"]->fill(deta_max);

    _h["m_emubblead_emu3b"]->fill(mass_emubb/GeV);

    _h["m_emubbfromtop_emu3b"]->fill(mass_emubbfromtop/GeV);

    _h["dR_emubblead_thirdbjet_emu3b"]->fill(dR_emubblead_thirdbjet);
    _h["dR_emubbfromtop_leadaddb_emu3b"]->fill(dR_emubbfromtop_leadadd);

    // light jets
    if (lightjets.size()>0) {
      _h["lead_lightjet_pt_emu3b"]->fill(lightjets[0].pT()/GeV);
      _h["lead_lightjet_eta_emu3b"]->fill(lightjets[0].abseta());

      const double dR_emubbfromtop_leadlightjet = deltaR(lljjfromtop, lightjets[0].mom());
      _h["dR_emubbfromtop_leadlightjet_emu3b"]->fill(dR_emubbfromtop_leadlightjet);
      _h["ptdiff_leadaddbjet_leadlightjet_emu3b"]->fill((lightjets[0].pT() - bfromQCD[0].pT())/GeV);
    }

    if (nbjets < 4) vetoEvent;

     //nlight jets
     _d["nlightjets_emu4b"]->fill(nlightjets>=3? "≥3"s : toString(nlightjets));
     // b-jet pTs
     _h["lead_bjet_pt_emu4b"]->fill(bjets[0].pT()/GeV);
     _h["sublead_bjet_pt_emu4b"]->fill(bjets[1].pT()/GeV);
     _h["third_bjet_pt_emu4b"]->fill(bjets[2].pT()/GeV);
     _h["fourth_bjet_pt_emu4b"]->fill(bjets[3].pT()/GeV);
     // bjets eta
     _h["lead_bjet_eta_emu4b"]->fill(bjets[0].abseta());
     _h["sublead_bjet_eta_emu4b"]->fill(bjets[1].abseta());
     _h["third_bjet_eta_emu4b"]->fill(bjets[2].abseta());
     _h["fourth_bjet_eta_emu4b"]->fill(bjets[3].abseta());

     // bjets assigned to top
     _h["lead_bjetfromtop_pt_emu4b"]->fill(bfromtop[0].pT()/GeV);
     _h["sublead_bjetfromtop_pt_emu4b"]->fill(bfromtop[1].pT()/GeV);
     _h["lead_bjetfromtop_eta_emu4b"]->fill(bfromtop[0].abseta());
     _h["sublead_bjetfromtop_eta_emu4b"]->fill(bfromtop[1].abseta());

     // additional b-jets
     _h["lead_bjetfromQCD_pt_emu4b"]->fill(bfromQCD[0].pT()/GeV);
     _h["sublead_bjetfromQCD_pt_emu4b"]->fill(bfromQCD[1].pT()/GeV);
     _h["lead_bjetfromQCD_eta_emu4b"]->fill(bfromQCD[0].abseta());
     _h["sublead_bjetfromQCD_eta_emu4b"]->fill(bfromQCD[1].abseta());

     // HT
     _h["ht_had_emu4b"]->fill(hthad/GeV);
     _h["ht_all_emu4b"]->fill(htall/GeV);

     // leading bb pair
     _h["m_bb_leading_emu4b"]->fill(jsum.mass()/GeV);
     _h["pt_bb_leading_emu4b"]->fill(jsum.pT()/GeV);
     _h["dR_bb_leading_emu4b"]->fill(dr_leading);

     // b-jets pair assigned to top
     _h["m_bbfromtop_emu4b"]->fill(jfromtopsum.mass()/GeV);
     _h["pt_bbfromtop_emu4b"]->fill(jfromtopsum.pT()/GeV);

     // closest bb pair
     _h["m_bb_closest_emu4b"]->fill(bb_closest.mass()/GeV);
     _h["pt_bb_closest_emu4b"]->fill(bb_closest.pT()/GeV);
     _h["dR_bb_closest_emu4b"]->fill(dr_closest);


     const FourMomentum bb_add = bfromQCD[0].mom() + bfromQCD[1].mom();
     _h["m_bbadd_emu4b"]->fill(bb_add.mass()/GeV);
     _h["pt_bbadd_emu4b"]->fill(bb_add.pT()/GeV);

     // global
     _h["dR_avg_emu4b"]->fill(dr_avg);
     _h["deta_jj_max_emu4b"]->fill(deta_max);
     _h["m_emubblead_emu4b"]->fill(mass_emubb/GeV);
     _h["m_emubbfromtop_emu4b"]->fill(mass_emubbfromtop/GeV);
     _h["dR_emubblead_thirdbjet_emu4b"]->fill(dR_emubblead_thirdbjet);
     _h["dR_emubbfromtop_leadaddb_emu4b"]->fill(dR_emubbfromtop_leadadd);

     if (lightjets.size()>0) {
       _h["lead_lightjet_pt_emu4b"]->fill(lightjets[0].pT()/GeV);
       _h["lead_lightjet_eta_emu4b"]->fill(lightjets[0].abseta());

       const double dR_emubbfromtop_leadlightjet = deltaR(lljjfromtop, lightjets[0].mom());
       _h["dR_emubbfromtop_leadlightjet_emu4b"]->fill(dR_emubbfromtop_leadlightjet);
       _h["ptdiff_leadaddbjet_leadlightjet_emu4b"]->fill((lightjets[0].pT() - bfromQCD[0].pT())/GeV);
     }

  }

  void finalize() {
    // Normalise all histograms
    normalize(_h, 1000.);
    normalize(_d, 1000.);
  }

  double drlb_min(const DressedLepton& lepton, const Jets& bjets) const {
    double mindr = 1e3;
    for (size_t i = 0; i < bjets.size(); ++i) {
      const double dr = deltaR(lepton, bjets[i]);
      if (dr < mindr)  mindr = dr;
    }
    return mindr;
  }

  /// @todo Replace with std::next_permutation when adopting C++20
  bool next_permutation(vector<Jet>::iterator first, vector<Jet>::iterator last) const {
    auto r_first = std::make_reverse_iterator(last);
    auto r_last = std::make_reverse_iterator(first);
    auto left = std::is_sorted_until(r_first, r_last,
                                     [](const Jet& l, const Jet& r) { return l.pT() > r.pT(); });

    if (left != r_last) {
      auto right = std::upper_bound(r_first, left, *left,
                                    [](const Jet& l, const Jet& r) { return l.pT() > r.pT(); });
      std::iter_swap(left, right);
    }

    std::reverse(left.base(), last);
    return left != r_last;
  }

  double bjetOriginWeight(const DressedLeptons& lepton, const Jets& bjets,
                          double mindRbb, double maxdRbb,
                          double mindRl1b, double mindRl2b) const {

    if (bjets.size() == 2) return 1.; //Return 1 if only two bjets are available

    if (lepton.size() < 2) return 0.;
    if (bjets.size() < 2)  return 0.;

    double diff_mindRl1b_dRl1b = mindRl1b - deltaR(lepton[0], bjets[0]);
    double diff_mindRl2b_dRl2b = mindRl2b - deltaR(lepton[1], bjets[1]);

    double weight = 1.0;
    if (bjets.size() < 4) {
      //dR between top and gluon b-jet
      const double dR_b1b3 = deltaR(bjets[0], bjets[2]);
      const double dR_b2b3 = deltaR(bjets[1], bjets[2]);
      const double maxdR_btop_bgluon = max(dR_b1b3, dR_b2b3);
      const double diff_maxdRbb_maxdRbb3 = maxdRbb - maxdR_btop_bgluon;
      weight *= exp(-1.*diff_mindRl1b_dRl1b*diff_mindRl1b_dRl1b);
      weight *= exp(-1.*diff_mindRl2b_dRl2b*diff_mindRl2b_dRl2b);
      weight *= exp(-1.*diff_maxdRbb_maxdRbb3*diff_maxdRbb_maxdRbb3);
    }
    else {
      const double dR_b3b4 = deltaR(bjets[2], bjets[3]);
      const double diff_mindRbb_dRb3b4 = mindRbb - dR_b3b4;
      weight *= exp(-1.*diff_mindRbb_dRb3b4*diff_mindRbb_dRb3b4);
      weight *= exp(-1.*diff_mindRl1b_dRl1b*diff_mindRl1b_dRl1b);
      weight *= exp(-1.*diff_mindRl2b_dRl2b*diff_mindRl2b_dRl2b);
    }

    return weight;
  }

 private:

  map<string,Histo1DPtr> _h;
  map<string,BinnedHistoPtr<string>> _d;

};

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2024_I2809112);

}
