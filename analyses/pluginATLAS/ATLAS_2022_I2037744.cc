// -*- C++ -*
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// Semileptonic ttbar single- and double-differential cross-sections with high pT top at 13 TeV
  class ATLAS_2022_I2037744 : public Analysis {
    public:

        /// Constructor
        RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2022_I2037744);

        void init() {

            // Define cut objects on eta, eta neutrino and leptons
            Cut eta_full =  Cuts::abseta < 5.0;
            Cut lep_cuts = Cuts::abseta < 2.5 && Cuts::pT > 27*GeV;

            // All final state particles
            const FinalState fs(eta_full);

            // Final state photons for loose lepton dressing for inputs to jets and MET
            IdentifiedFinalState all_photons(fs, PID::PHOTON);

            // Final state photons, not from taus, for analysis lepton dressing
            PromptFinalState photons(all_photons, TauDecaysAs::NONPROMPT);
            declare(photons, "photons");

            // Final state electrons, including from prompt tau decays
            PromptFinalState electrons(Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);
            declare(electrons, "electrons");

            // Analysis dressed electrons
            LeptonFinder dressedelectrons(electrons, photons, 0.1, lep_cuts);
            declare(dressedelectrons, "dressedelectrons");

            // "All" dressed electrons to be removed from input to jetbuilder
            LeptonFinder ewdressedelectrons(electrons, all_photons, 0.1, eta_full);
            declare(ewdressedelectrons, "ewdressedelectrons");

            //Final state muons, including from prompt tau decays
            PromptFinalState muons(Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT);
            declare(muons, "muons");

            //Analysis dressed muons
            LeptonFinder dressedmuons(muons, photons, 0.1, lep_cuts);
            declare(dressedmuons, "dressedmuons");

            //"All" dressed muons to be removed from input to jetbuilder and for use in METbuilder
            LeptonFinder ewdressedmuons(muons, all_photons, 0.1, eta_full);
            declare(ewdressedmuons, "ewdressedmuons");

            //Neutrinos to be removed from input to jetbuilder, acceptTauDecays=true
            IdentifiedFinalState nu_id;
            nu_id.acceptNeutrinos();
            PromptFinalState neutrinos(nu_id, TauDecaysAs::PROMPT);

            //Small-R jets
            VetoedFinalState vfs(fs);
            vfs.addVetoOnThisFinalState(ewdressedelectrons);
            vfs.addVetoOnThisFinalState(ewdressedmuons);
            vfs.addVetoOnThisFinalState(neutrinos);
            FastJets jets(vfs, JetAlg::ANTIKT, 0.4);
            jets.useInvisibles(JetInvisibles::DECAY);
            declare(jets, "jets");

            //MET
            declare(MissingMomentum(), "MissingMomentum");

            // External bins for 2D plots
            vector<double> n_jet_2D_bins = {0.5,1.5,2.5,10.0}; //Extra wide final bin mistakenly used in analysis, replicated here
            vector<double> Top_boosted_rc_pt_2D_bins = {355.0,398.0,496.0,2000.0};

            //Book Histograms using custom function (handles HEPData offset and relative hists)
            book_hist("sigma_ttbar",1,false);
            book_hist("Top_boosted_rc_pt",2);
            book_hist("Top_boosted_leptonic_pt",5);
            book_hist("ttbar_boosted_rc_m",8);
            book_hist("hadTop_boosted_rc_y",11);
            book_hist("lepTop_boosted_y",14);
            book_hist("ttbar_boosted_rc_y",17);
            book_hist("boosted_rc_HT",20);
            book_hist("dphi_lepb_hadTop",23);
            book_hist("ttbar_boosted_rc_pt",26);
            book_hist("dphi_hadTop_lepTop",29);
            book_hist("HTall",32);
            book(_njets, 36,1,1);
            book_hist("LeadAddJet_pt",38);
            book_hist("LeadAddJet_hadTop_m",41);
            book_hist("dphi_LeadAddJet_hadTop",44);
            book_hist("dphi_SubLeadAddJet_hadTop",47);
            book_hist("dphi_LeadAddJet_SubLeadAddJet",50);
            book_hist("SubLeadAddJet_pt",53);

            book_2Dhist("LeadAddJet_pt_2D_Nextrajets",n_jet_2D_bins,56);
            book_2Dhist("LeadAddJet_pt_2D_Top_boosted_rc_pt",Top_boosted_rc_pt_2D_bins,68);
            book_2Dhist("dphi_LeadAddJet_hadTop_2D_Top_boosted_rc_pt",Top_boosted_rc_pt_2D_bins,80);
            book_2Dhist("dphi_LeadAddJet_hadTop_2D_Nextrajets",n_jet_2D_bins,92);
        }

        void analyze(const Event& event) {

            //----------Projections
            DressedLeptons electrons = apply<LeptonFinder>(event, "dressedelectrons").dressedLeptons();
            DressedLeptons muons = apply<LeptonFinder>(event, "dressedmuons").dressedLeptons();

            // We a need seperate jet collection with 25GeV cut to perform OR with leptons
            const Jets& smalljets_25 = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25.0*GeV && Cuts::abseta <= 2.5);

            idiscardIfAnyDeltaRLess(electrons, smalljets_25, 0.4);
            idiscardIfAnyDeltaRLess(muons, smalljets_25, 0.4);

            // Analysis small-R jets, 26GeV cut
            const Jets& smalljets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 26.0*GeV && Cuts::abseta <= 2.5);
            // Need jets for re-clustering with a 30GeV cut (default in AnalysisTop), so make a seperate collection
            Jets smalljets_for_rc;

            map<double,bool> b_tagging_info_small_jets; // Need to append b-tagging information to pseudojets for top-tagging to work, use momentum<->btag map
            map<double,bool> b_tagging_info_rc_jets;
            map<double,bool> Addjet_veto_info; // Also append information about which small-R jets are subjets of the chosen top-quark or the leptonic b-jet, again use momentum<->subjet map

            for (Jet jet : smalljets) {
	      b_tagging_info_small_jets[jet.pt()]=jet.bTagged(Cuts::pT > 5.0*GeV);
	      if (jet.pt() >= 30.0*GeV){
		smalljets_for_rc += jet;
	      }
            }

            const FourMomentum met = apply<MissingMomentum>(event, "MissingMomentum").missingMomentum();

            // Define reclustered jets AntiKT 1.0
            fastjet::ClusterSequence antikt_cs(smalljets_for_rc, fastjet::JetDefinition(fastjet::antikt_algorithm, 1.0));
            PseudoJets reclustered_jets = antikt_cs.inclusive_jets();

            // trim the jets
            PseudoJets TrimmedReclusteredJets;

	    /// @todo Store rather than rebuild for every event
            fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, 0.01), fastjet::SelectorPtFractionMin(0.05));
            for (PseudoJet pjet : reclustered_jets) {
	      fastjet::PseudoJet candidate_trim = trimmer(pjet);
	      bool b_tagged = false;
	      std::vector<fastjet::PseudoJet> constituents = candidate_trim.constituents();
	      for (unsigned int iCons = 0; iCons<constituents.size(); iCons++) {
		if (b_tagging_info_small_jets[constituents[iCons].pt()]) b_tagged = true;
	      }
	      FourMomentum trfj_mom = momentum(candidate_trim);
	      if (trfj_mom.pt() <= 355*GeV)  continue;
	      if (trfj_mom.abseta() < 2.0) {
		TrimmedReclusteredJets.push_back(candidate_trim);
		b_tagging_info_rc_jets[candidate_trim.perp()]=b_tagged;
	      }
            }
            TrimmedReclusteredJets = fastjet::sorted_by_pt(TrimmedReclusteredJets);
	    


            //----------Event selection

            // SINGLE LEPTON
            bool single_electron=(electrons.size() == 1) && (muons.empty());
            bool single_muon=(muons.size() == 1) && (electrons.empty());
            DressedLepton* lepton = NULL;
            if (single_electron) {
                lepton = &electrons[0];
            }
            else if (single_muon) {
                lepton = &muons[0];
            }
            else {
                vetoEvent;
            }

            //MET
            if (met.pt()<20.0*GeV) vetoEvent;
            double transmass = mT(momentum(*lepton), met);
            //MET+MWT
            if ((met.pt()+transmass)<60.0*GeV) vetoEvent;

            //SMALL-R JET MULTIPLICITY
            if (smalljets.size()<2) vetoEvent;
            if (TrimmedReclusteredJets.size()==0) vetoEvent;

            //TOP-TAGGED RC JET
            PseudoJet HadTopJet;
            bool ThereIsHadTop = false;
            for (PseudoJet rc_jet : TrimmedReclusteredJets){
                FourMomentum rc_jet_mom = momentum(rc_jet);
                double dR_lepJet = deltaR(rc_jet_mom,momentum(*lepton));
                if (single_electron && dR_lepJet < 1.) continue;
                if (!b_tagging_info_rc_jets[rc_jet_mom.pt()]) continue;
                if (rc_jet_mom.mass() > 120.0*GeV && rc_jet_mom.mass() < 220.0*GeV) {
                    HadTopJet = rc_jet;
                    ThereIsHadTop = true;
                    break; //Pick highest pT trimmed RC jet passing requirements
                }
            }
            if (!ThereIsHadTop) vetoEvent;


            // BTAGGED JET ON LEPTONIC SIDE
            Jet LepbJet;
            double smallest_dR_bjetlep=2.0;
            bool ThereIsLepbJet = false;
            for (Jet jet : smalljets) {
                // leptonic bjet cannot be constituent of top-jet
                std::vector<fastjet::PseudoJet> constituents = HadTopJet.constituents();
                bool issubjet=false;
                for(PseudoJet subjet : constituents) {
                    // can't do direct equality because smalljets and RCsubjets are different jet collections, so we do an ugly pT comparison
                    if (jet.pt() == momentum(subjet).pt() ) {
                        issubjet=true;
                        Addjet_veto_info[jet.pt()] = true;
                    }
                }
                if (issubjet) continue;
                if (!b_tagging_info_small_jets[jet.pt()]) continue; // Must be b-tagged (do after so we can also fill addjet veto in same loop)

                const double dR_bjetlep = deltaR(jet, *lepton);
                if (dR_bjetlep > smallest_dR_bjetlep) continue;
                else {
                    smallest_dR_bjetlep = dR_bjetlep;
                    LepbJet = jet; //Take b-tagged non-subjet small-R jet closest in dR to lepton
                    ThereIsLepbJet = true;
                }
            }
            if (!ThereIsLepbJet) vetoEvent;
            Addjet_veto_info[momentum(LepbJet).pt()] = true;

            // MLB
            double mlb = (lepton->momentum() + LepbJet.momentum()).mass();
            if (mlb >= 180.0*GeV) vetoEvent;

            // Reconstruct leptonically decaying top-jet
            const double nu_pz = computeneutrinoz(lepton->momentum(), met, LepbJet.momentum());
            FourMomentum neutrino( sqrt(sqr(met.px()) + sqr(met.py()) + sqr(nu_pz)), met.px(), met.py(), nu_pz);
            FourMomentum LeptonicTop = lepton->momentum() + neutrino + LepbJet.momentum();
            FourMomentum HadronicTop = momentum(HadTopJet);
            FourMomentum pttbar = HadronicTop + LeptonicTop;

            // Lastly find additional jets
            double HT_all = 0.0; //Set up variable for pT sum of ttbar and all additional jets
            Jets addJets;
            for (const Jet& jet : smalljets) {
                // ignore all sub-jets of hadronic top and the b-tagged jet on the leptonic side
                if (Addjet_veto_info[jet.pt()]) continue;
                addJets += jet;
                HT_all = HT_all + jet.pt();
            }
            FourMomentum leading_addjet;
            FourMomentum subleading_addjet;
            FourMomentum p_hadtop_leading_addjet;// = HadronicTop;
            if (addJets.size() > 0) {
                leading_addjet = addJets[0].momentum();
                p_hadtop_leading_addjet = HadronicTop + leading_addjet;
                if (addJets.size() > 1) {
                    subleading_addjet = addJets[1].momentum();
                }
            }

            // calculate some observables
            const double HT_ttbar = HadronicTop.pt() + LeptonicTop.pt();
            HT_all = HT_all + HT_ttbar;
            const double dphi_lepb_hadTop = deltaPhi(LepbJet.momentum(), HadronicTop)/PI;
            const double dphi_hadTop_lepTop = deltaPhi(HadronicTop, LeptonicTop)/PI;

            // Observables
            fillHist("sigma_ttbar",             0.5,false);
            fillHist("Top_boosted_rc_pt",       HadronicTop.pt()/GeV);
            fillHist("Top_boosted_leptonic_pt", LeptonicTop.pt()/GeV);
            fillHist("ttbar_boosted_rc_m",      pttbar.mass()/GeV);
            fillHist("hadTop_boosted_rc_y",     HadronicTop.absrap());
            fillHist("lepTop_boosted_y",        LeptonicTop.absrap());
            fillHist("ttbar_boosted_rc_y",      pttbar.absrap());
            fillHist("boosted_rc_HT",           HT_ttbar/GeV);
            fillHist("dphi_lepb_hadTop",        dphi_lepb_hadTop);
            fillHist("ttbar_boosted_rc_pt",     pttbar.pt()/GeV);
            fillHist("dphi_hadTop_lepTop",      dphi_hadTop_lepTop);
            fillHist("HTall",                   HT_all/GeV);
            _njets->fill(map2string(min(addJets.size(), 6u)));

            if (addJets.size() > 0) {
                const double dphi_leadaddjet_hadTop = deltaPhi( leading_addjet,HadronicTop ) / PI;
                fillHist("LeadAddJet_pt",          leading_addjet.pt()/GeV);
                fillHist("LeadAddJet_hadTop_m",    p_hadtop_leading_addjet.mass()/GeV);
                fillHist("dphi_LeadAddJet_hadTop", dphi_leadaddjet_hadTop);

                // 2D Observables
                fillHist2D("LeadAddJet_pt_2D_Nextrajets", min(addJets.size(), 6u), leading_addjet.pt()/GeV);
                fillHist2D("LeadAddJet_pt_2D_Top_boosted_rc_pt", HadronicTop.pt()/GeV, leading_addjet.pt()/GeV);
                fillHist2D("dphi_LeadAddJet_hadTop_2D_Top_boosted_rc_pt", HadronicTop.pt()/GeV, dphi_leadaddjet_hadTop);
                fillHist2D("dphi_LeadAddJet_hadTop_2D_Nextrajets", min(addJets.size(), 6u), dphi_leadaddjet_hadTop);
            }

            if (addJets.size() > 1) {
                const double dphi_subleadaddjet_hadTop = deltaPhi( subleading_addjet, HadronicTop ) / PI;
                const double dphi_leadaddjet_subleadaddjet = deltaPhi( leading_addjet, subleading_addjet ) / PI;
                fillHist("dphi_SubLeadAddJet_hadTop",     dphi_subleadaddjet_hadTop );
                fillHist("dphi_LeadAddJet_SubLeadAddJet", dphi_leadaddjet_subleadaddjet );
                fillHist("SubLeadAddJet_pt",              subleading_addjet.pt()/GeV);
            }

        }

      
        void finalize() {

            // Normalize to cross-section
            const double sf = crossSection() / picobarn / sumOfWeights();

            for (auto& hist : _h) {
              scale(hist.second, sf);
              if (hist.first.find("_norm") != string::npos)  normalize(hist.second, 1.0, false);
            }
            scale(_njets, sf);
            for (auto& hist : _h_multi) {
              scale(hist.second, sf);
              if (hist.first.find("_norm") != string::npos) {
                normalize(hist.second, 1.0, false);
              }
              divByGroupWidth(hist.second);
            }

        }

      
    private:

        // HepData entry has dummy "Table of Contents", for both 1D and 2D hists need to offset tables by one unit
        void book_hist(const string& name, unsigned int table, bool do_norm = true) {
          book(_h[name], table+1, 1, 1);
          if (do_norm) {
            book(_h[name+"_norm"], table+3, 1, 1);
          }
        }

        void book_2Dhist(const string& name, const std::vector<double>& doubleDiff_bins, unsigned int table) {
          book(_h_multi[name+"_norm"], doubleDiff_bins);
          book(_h_multi[name], doubleDiff_bins);
          for (size_t i=0; i < _h_multi[name]->numBins(); ++i) {
            book(_h_multi[name+"_norm"]->bin(i+1), table+i+1, 1, 1);
            book(_h_multi[name]->bin(i+1), table+i+4, 1, 1);
          }
        }

        // Fill abs and nomralised hists at same time
        void fillHist(const string& name, double value, bool do_norm = true) {
          _h[name]->fill(value);
          if (do_norm)  _h[name+"_norm"]->fill(value);
        }

        void fillHist2D(const string& name, double externalbin, double val) {
          _h_multi[name]->fill(externalbin, val);
          _h_multi[name+"_norm"]->fill(externalbin, val);
        }

        string map2string(const size_t njets) const {
          if (njets == 0)  return "0";
          if (njets == 1)  return "1";
          if (njets == 2)  return "2";
          if (njets  < 5)  return "3.0 - 4.0";
          return "$>$4";
        }


        double computeneutrinoz(const FourMomentum& lepton, const FourMomentum& met, const FourMomentum& lbjet) const {
            const double m_W = 80.385*GeV; // mW
            const double beta = m_W*m_W - lepton.mass()*lepton.mass() + 2.0*lepton.px()*met.px() + 2.0*lepton.py()*met.py();
            const double delta = lepton.E()*lepton.E()*( beta*beta + (2.0*lepton.pz()*met.pt())*(2.0*lepton.pz()*met.pt()) - (2.0*lepton.E()*met.pt())*(2.0*lepton.E()*met.pt()) );
            if(delta <= 0) {
                //imaginary solution, use real part
                double pzneutrino = 0.5*lepton.pz()*beta / (lepton.E()*lepton.E() - lepton.pz()*lepton.pz());
                return pzneutrino;
            }
            double pzneutrinos[2] = {0.5 * (lepton.pz()*beta + sqrt(delta)) / (lepton.E()*lepton.E() - lepton.pz()*lepton.pz()),
                                     0.5 * (lepton.pz()*beta - sqrt(delta)) / (lepton.E()*lepton.E() - lepton.pz()*lepton.pz())};
            FourMomentum topCands[2];
            for(int i=0; i<2; ++i) {
                FourMomentum neutrino;
                neutrino.setXYZM( met.px(), met.py(), pzneutrinos[i], 0.0 );
                topCands[i] = neutrino + lbjet + lepton;
            }
            // Pick neutrino solution that results in smallest top mass
            if (topCands[0].mass() <= topCands[1].mass() ) {
                return pzneutrinos[0];
            } else {
                return pzneutrinos[1];
            }
        }

        // Histogram pointer maps
        map<string, Histo1DPtr> _h;
        BinnedHistoPtr<string> _njets;
        map<string, Histo1DGroupPtr> _h_multi;
  };


  RIVET_DECLARE_PLUGIN(ATLAS_2022_I2037744);

}
