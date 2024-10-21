// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief Z+b(b) and Z+c at 13 TeV
  class ATLAS_2024_I2771257 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2024_I2771257);

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      _mode = 0;
      if ( getOption("LMODE") == "EL" ) _mode = 1;
      if ( getOption("LMODE") == "MU" ) _mode = 2;

      const FinalState fs;
      // Define fiducial cuts for the leptons in the ZFinder
      const Cut lepcuts = (Cuts::pT > 27*GeV) && (Cuts::abseta < 2.5);
      DileptonFinder zfinderE(91.2*GeV, 0.1, lepcuts && Cuts::abspid == PID::ELECTRON, Cuts::massIn(76*GeV, 106*GeV));
      DileptonFinder zfinderM(91.2*GeV, 0.1, lepcuts && Cuts::abspid == PID::MUON, Cuts::massIn(76*GeV, 106*GeV));
      declare(zfinderE, "zfinderE");
      declare(zfinderM, "zfinderM");

      //Build jets using ATLAS AntiKt4WZtruthJets definition removing dressed W, Z leptons

      // Photons
      FinalState photons(Cuts::abspid == PID::PHOTON);

      // Muons
      PromptFinalState bare_mu(Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT);
      LeptonFinder all_dressed_mu(bare_mu, photons, 0.1, Cuts::abseta < 2.5);

      // Electrons
      PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);
      LeptonFinder all_dressed_el(bare_el, photons, 0.1, Cuts::abseta < 2.5);

      //Jet forming
      VetoedFinalState vfs(FinalState(Cuts::abseta < 4.5));
      vfs.addVetoOnThisFinalState(all_dressed_el);
      vfs.addVetoOnThisFinalState(all_dressed_mu);

      FastJets jetfs(vfs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::DECAY);
      declare(jetfs, "jets");

      //HF hadrons for b and c jet def
      declare(HeavyHadrons(), "HFHadrons");

      // fiducial1B
      book(_h["fiducial1B_Z_Pt"],         4, 1, 1);
      book(_h["fiducial1B_leadBJet_Pt"],  5, 1, 1);
      book(_h["fiducial1B_ZleadBJet_DR"], 6, 1, 1);

      // fiducial2B
      book(_h["fiducial2B_diBJets_DPhi"], 7, 1, 1);
      book(_h["fiducial2B_diBJets_M"],    8, 1, 1);

      // fiducial1C
      book(_h["fiducial1C_Z_Pt"],         9, 1, 1);
      book(_h["fiducial1C_leadCJet_Pt"], 10, 1, 1);
      book(_h["fiducial1C_leadCJet_xF"], 11, 1, 1);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // --------------------------------------- object sel: leptons -----------------------------------------------

      const DileptonFinder& zfinderE = apply<DileptonFinder>(event, "zfinderE");
      const Particles& els = zfinderE.constituents();
      const DileptonFinder& zfinderM = apply<DileptonFinder>(event, "zfinderM");
      const Particles& mus = zfinderM.constituents();

     // default is to run average of Z->ee and Z->mm
     // use LMODE option to pick one channel
      if ( (els.size() + mus.size()) != 2 )  vetoEvent;

      if (      _mode == 0 && !(els.size()==2 || mus.size()==2) )  vetoEvent;
      else if ( _mode == 1 && !(els.size() == 2 && mus.empty()) )  vetoEvent;
      else if ( _mode == 2 && !(els.empty() && mus.size() == 2) )  vetoEvent;

      double Z_pT = 0., Zrap = 0., Zphi = 0.;
      if ( els.size()==2 ) {
        Z_pT = zfinderE.boson().pt()/GeV;
        Zphi = zfinderE.boson().phi();
        Zrap = zfinderE.boson().rapidity();
      } else {
        Z_pT = zfinderM.boson().pt()/GeV;
        Zphi = zfinderM.boson().phi();
        Zrap = zfinderM.boson().rapidity();
      }

      // --------------------------------------- object sel: jets ---------------------------------------------------
      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::absrap < 2.5);

      // --------------------------------------- object sel: OR ---------------------------------------------------
      //OR 0.4 between jets and dressed leptons
      idiscardIfAnyDeltaRLess(jets, els, 0.4);
      idiscardIfAnyDeltaRLess(jets, mus, 0.4);

      //b- and c- jet definition according to cone R=0.4 labelling
      Jets bjets, cjets;
      const HeavyHadrons &ha = apply<HeavyHadrons>(event,"HFHadrons");
      const Particles allBs = ha.bHadrons(Cuts::pT > 5*GeV);
      const Particles allCs = ha.cHadrons(Cuts::pT > 5*GeV);
      Particles matchedBs, matchedCs;

      std::vector<int> bJet_idx_v;
      int jet_idx =0;

      for (const Jet& j : jets) {
        Jet closest_j;
        Particle closest_b;
        double minDR_j_b = 10;

        for (const Particle& bHad : allBs) {
          bool alreadyMatched = false;
          for (const Particle& bMatched : matchedBs) {
            alreadyMatched |= bMatched.isSame(bHad);
          }
          if (alreadyMatched)  continue;

          double DR_j_b = deltaR(j, bHad);
          if (DR_j_b <= 0.3 && DR_j_b < minDR_j_b) {
            minDR_j_b = DR_j_b;
            closest_j = j;
            closest_b = bHad;
          }
        }

        if (minDR_j_b < 0.3) {
          bjets += closest_j;
          matchedBs += closest_b;
          bJet_idx_v.push_back(jet_idx);
        }
        ++jet_idx;
      }

      jet_idx = 0;
      for (const Jet& j : jets) {
        Jet closest_j;
        Particle closest_c;
        double minDR_j_c = 10;

        for (const Particle& cHad : allCs) {
          bool alreadyMatched = false;
          for (const Particle& cMatched : matchedCs) {
            alreadyMatched |= cMatched.isSame(cHad);
          }

          // Do not match c-jets if already b-matched
          for (size_t bJet_i  = 0; bJet_i<bJet_idx_v.size(); ++bJet_i) {
            if (jet_idx == bJet_idx_v[bJet_i]) alreadyMatched = true;
          }

          if (alreadyMatched)  continue;

          double DR_j_c = deltaR(j, cHad);
          if ( DR_j_c <= 0.3 && DR_j_c < minDR_j_c) {
            minDR_j_c = DR_j_c;
            closest_j = j;
            closest_c = cHad;
          }
        }

        if (minDR_j_c < 0.3) {
          cjets += closest_j;
          matchedCs += closest_c;
        }
        ++jet_idx;
      }

      // --------------------------------------- event sel: Nbjets ------------------------------------------------------
      if (bjets.size() < 1 && cjets.size() < 1)  vetoEvent;

      // --------------------------------------- event classification ------------------------------------------------------
      // string region = "";
      bool pass_fid1B = false;
      bool pass_fid1C = false;
      bool pass_fid2B = false;

      // based on the flavour of the leading HF jet
      if ( (bjets.size()>=1 && cjets.size()==0) ||
           (bjets.size()>=1 && cjets.size()>=1 && bjets[0].pT() > cjets[0].pT()) ) {
        pass_fid1B = true;
      }
      if ( (cjets.size()>=1 && bjets.size()==0) ||
           (cjets.size()>=1 && bjets.size()>=1 && bjets[0].pT() < cjets[0].pT()) ) {
        pass_fid1C = true;
      }

      // based on the flavour of the leading HF jet
      if ( (bjets.size()>=2 && cjets.size()==0) ||
           (bjets.size()>=2 && cjets.size()>=1 && bjets[0].pT() > cjets[0].pT()) ) {
        pass_fid2B = true;
      }
      // -------------------------------------------------------------------------------------------------------------------



      // --------------------------------------- observable cal -------------------------------------------------------

      double leadBJet_Pt = -100, ZleadBJet_DR = -100, dRapZb = -100, dPhiZb = -100, diBJets_M = -100, diBJets_DPhi = -100;

      if (pass_fid1B){
        leadBJet_Pt = bjets[0].pT();
        dRapZb=fabs(Zrap-bjets[0].rap());
        dPhiZb=acos(cos(fabs(Zphi-bjets[0].phi())));
        ZleadBJet_DR = sqrt(dRapZb*dRapZb+dPhiZb*dPhiZb);
      }

      if (pass_fid2B){
        diBJets_M = (bjets[0].momentum() + bjets[1].momentum()).mass();
        diBJets_DPhi = acos(cos(fabs(bjets[0].phi()-bjets[1].phi())));
      }

      // Feynam scaling variable defined as (2 pT*sinh(eta))/sqrt(s) will have the negative range
      // and pT*sinh(eta) == pZ, which is not exactly accurate for massive objects, such as jet
      double leadCJet_Pt = -100, leadCJet_xF = -100;

      if (pass_fid1C){
        leadCJet_Pt = cjets[0].pT();
        leadCJet_xF = 2*fabs(cjets[0].pz()/GeV)/13000.;
      }

      // ------------------------------------------ hist fill ---------------------------------------------------------
      if (pass_fid1B){
        _h["fiducial1B_Z_Pt"]->fill(Z_pT/GeV);
        _h["fiducial1B_leadBJet_Pt"]->fill(leadBJet_Pt/GeV);
        _h["fiducial1B_ZleadBJet_DR"]->fill(ZleadBJet_DR);
      }

      if (pass_fid2B){
        _h["fiducial2B_diBJets_M"]->fill(diBJets_M/GeV);
        _h["fiducial2B_diBJets_DPhi"]->fill(diBJets_DPhi);
      }

      if (pass_fid1C){
        _h["fiducial1C_Z_Pt"]->fill(Z_pT/GeV);
        _h["fiducial1C_leadCJet_Pt"]->fill(leadCJet_Pt/GeV);
        _h["fiducial1C_leadCJet_xF"]->fill(leadCJet_xF);
      }

      // -------------------------------------------------------------------------------------------------------------

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = _mode? 1.0 : 0.5;
      scale(_h, sf * crossSectionPerEvent());
    }

  private:

     size_t _mode;

    map<string, Histo1DPtr> _h;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2024_I2771257);

}
