// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Math/LorentzTrans.hh"

namespace Rivet {


  /// @brief All-hadronic ttbar at 13 TeV
  class ATLAS_2022_I2077575 : public Analysis {
  public:

      /// Constructor
      RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2022_I2077575);

      /// Book histograms and initialise projections before the run
      void init() {

        // Get options particle-level only.
        _mode = 0;
        if ( getOption("TMODE") == "PARTICLE" ) _mode = 0;
        if ( getOption("TMODE") == "BOTH" ) _mode = 1;

        // External bins for 2D and 3D cross-sections
        std::vector<double> t1_pt_2D_bins_1 = {0.5, 0.55, 0.6, 0.75, 2.0};
        std::vector<double> t1_pt_2D_bins_2 = {0.5, 0.55, 0.625, 0.75, 2.0};
        std::vector<double> t_and_tt_y_2D_bins = {0.0, 0.2, 0.5, 1.0, 2.0};
        std::vector<double> tt_pt_2D_bins = {0.0, 0.1, 0.2, 0.35, 1.0};
        std::vector<double> tt_m_3D_bins = {0.9, 1.2, 1.5, 4.0};

        //histogram booking
        book(_h["inclusive_particle"], 2, 1, 1);
        if (_mode)  book(_h["inclusive_parton"], 147, 1, 1);
        book_hist("t_pt", 	       3);
        book_hist("t_y",  	       4);
        book_hist("t1_pt",         5);
        book_hist("t1_y",          6);
        book_hist("t2_pt",         7);
        book_hist("t2_y",          8);
        book_hist("tt_m",          9);
        book_hist("tt_pt",         10);
        book_hist("tt_y",          11);
        book_hist("tt_chi",        12);
        book_hist("tt_yboost",     13);
        book_hist("tt_pout",       14);
        book_hist("tt_dPhi",       15);
        book_hist("tt_Ht",         16);
        book_hist("tt_cosThStar",  17);
        book_hist_2D("t1_pt_t2_pt_2D", 		    t1_pt_2D_bins_1, 	    18);
        book_hist_2D("t1_y_t2_y_2D", 		      t_and_tt_y_2D_bins, 	22);
        book_hist_2D("t1_y_t1_pt_2D", 		    t_and_tt_y_2D_bins, 	26);
        book_hist_2D("t2_y_t2_pt_2D", 		    t_and_tt_y_2D_bins, 	30);
        book_hist_2D("t1_pt_tt_pt_2D", 		    t1_pt_2D_bins_2, 	    34);
        book_hist_2D("t1_pt_tt_m_2D", 		    t1_pt_2D_bins_2, 	    38);
        book_hist_2D("tt_y_t1_pt_2D", 		    t_and_tt_y_2D_bins, 	42);
        book_hist_2D("tt_y_t1_y_2D", 		      t_and_tt_y_2D_bins, 	46);
        book_hist_2D("t1_y_tt_m_2D", 		      t_and_tt_y_2D_bins, 	50);
        book_hist_2D("tt_y_tt_m_2D", 		      t_and_tt_y_2D_bins, 	54);
        book_hist_2D("tt_pt_tt_m_2D", 		    tt_pt_2D_bins, 		    58);
        book_hist_2D("tt_y_tt_pt_2D",         t_and_tt_y_2D_bins, 	62);
        book_hist_2D("tt_y_1_tt_m_t1_pt_3D", 	tt_m_3D_bins, 		    66);
        book_hist_2D("tt_y_2_tt_m_t1_pt_3D", 	tt_m_3D_bins, 		    69);
        book_hist_2D("tt_y_3_tt_m_t1_pt_3D", 	tt_m_3D_bins, 		    72);

        // Projections
        const Cut dressed_lep = Cuts::abseta < 2.5 && Cuts::pT >= 25*GeV;
        const Cut all_dressed_lep = Cuts::abseta < 2.5;
        const Cut eta_full = Cuts::abseta < 4.5;

        // All final state particles
        const FinalState fs(eta_full);

        // Get photons to dress leptons
        const FinalState photons(Cuts::abspid == PID::PHOTON);

        // Projection to find the electrons
        PromptFinalState electrons(Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);
        LeptonFinder dressedelectrons(electrons, photons, 0.1, dressed_lep);
        declare(dressedelectrons, "elecs");
        LeptonFinder alldressedelectrons(electrons, photons, 0.1, all_dressed_lep);

        // Projection to find the muons
        PromptFinalState muons(Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT);
        LeptonFinder dressedmuons(muons, photons, 0.1, dressed_lep);
        declare(dressedmuons, "muons");
        LeptonFinder alldressedmuons(muons, photons, 0.1, all_dressed_lep);

        // Small-R jet clustering
        VetoedFinalState vfs(fs);
        vfs.addVetoOnThisFinalState(alldressedelectrons);
        vfs.addVetoOnThisFinalState(alldressedmuons);
        FastJets sjets(vfs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::DECAY);
        declare(sjets, "sjets");

        // Large-R jet clustering.
        FastJets ljets(fs, JetAlg::ANTIKT, 1.0, JetMuons::NONE, JetInvisibles::NONE);
        ljets.addTrf(new fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2),
                                         fastjet::SelectorPtFractionMin(0.05)));
        declare(ljets, "ljets");

        if (_mode) {
          PartonicTops partonTops;
          declare(partonTops, "partonicTops");
        }
      }


      void analyze(const Event& event) {

        if (_mode) {

          // Parton-level top quarks
          const Particles partonicTops = apply<PartonicTops>( event, "partonicTops").particlesByPt();
          FourMomentum top, tbar;
          bool foundT = false, foundTBar = false;
          for (const Particle& ptop : partonicTops) {
            const int pid = ptop.pid();
            if (pid == PID::TQUARK) {
              top = ptop.momentum();
              foundT = true;
            }
            else if (pid == -PID::TQUARK) {
              tbar = ptop.momentum();
              foundTBar = true;
            }
          }

          FourMomentum t1_parton, t2_parton, ttbar_parton;
          if ( foundT && foundTBar ) {
            t1_parton = top.pT2() > tbar.pT2() ? top : tbar;
            t2_parton = top.pT2() > tbar.pT2() ? tbar : top;
            ttbar_parton = t1_parton + t2_parton;

            if ( t1_parton.pT() > 500*GeV && t2_parton.pT() > 350*GeV) {

              const double chi_parton = calcChi(t1_parton, t2_parton);
              const double cosThetaStar_parton = abs(calcCosThetaStar(t1_parton, t2_parton));
              if (cosThetaStar_parton == -99) {
                MSG_DEBUG("ttbar going faster than light! Vetoing event. Try turning of partonic tops?");
                vetoEvent;
              }
              const double pout_parton = abs(calcPout(t1_parton, t2_parton));
              const double dPhi_parton = deltaPhi(t1_parton, t2_parton);

              const FourMomentum& randomTopParton = t1_parton; // : t2_parton;

              if (_mode) _h["inclusive_parton"]->fill(0);

              fill_hist_parton("t_pt", randomTopParton.pT()/TeV);
              fill_hist_parton("t_y",  randomTopParton.absrap());

              fill_hist_parton("t1_pt", t1_parton.pT()/TeV);
              fill_hist_parton("t1_y",  t1_parton.absrap());
              fill_hist_parton("t2_pt", t2_parton.pT()/TeV);
              fill_hist_parton("t2_y",  t2_parton.absrap());

              fill_hist_parton("tt_m",  ttbar_parton.mass()/TeV);
              fill_hist_parton("tt_pt", ttbar_parton.pT()/TeV);
              fill_hist_parton("tt_Ht", (t1_parton.pT() + t2_parton.pT())/TeV);
              fill_hist_parton("tt_y",  ttbar_parton.absrap());

              fill_hist_parton("tt_yboost", 0.5 * abs(t1_parton.rapidity() + t2_parton.rapidity()));
              fill_hist_parton("tt_chi", chi_parton);
              fill_hist_parton("tt_cosThStar", cosThetaStar_parton);
              fill_hist_parton("tt_pout", pout_parton/TeV);
              fill_hist_parton("tt_dPhi", dPhi_parton);

              fill_hist_2D_parton("t1_pt_t2_pt_2D", t1_parton.pT()/TeV, t2_parton.pT()/TeV);
              fill_hist_2D_parton("t1_y_t2_y_2D", t1_parton.absrap(), t2_parton.absrap());
              fill_hist_2D_parton("t1_y_t1_pt_2D", t1_parton.absrap(), t1_parton.pT()/TeV);
              fill_hist_2D_parton("t2_y_t2_pt_2D", t2_parton.absrap(), t2_parton.pT()/TeV);
              fill_hist_2D_parton("t1_pt_tt_pt_2D", t1_parton.pT()/TeV, ttbar_parton.pT()/TeV);
              fill_hist_2D_parton("t1_pt_tt_m_2D", t1_parton.pT()/TeV, ttbar_parton.mass()/TeV);
              fill_hist_2D_parton("tt_y_t1_pt_2D", ttbar_parton.absrap(), t1_parton.pT()/TeV);
              fill_hist_2D_parton("tt_y_t1_y_2D", ttbar_parton.absrap(), t1_parton.absrap());
              fill_hist_2D_parton("t1_y_tt_m_2D", t1_parton.absrap(), ttbar_parton.mass()/TeV);
              fill_hist_2D_parton("tt_y_tt_m_2D", ttbar_parton.absrap(), ttbar_parton.mass()/TeV);
              fill_hist_2D_parton("tt_pt_tt_m_2D", ttbar_parton.pT()/TeV, ttbar_parton.mass()/TeV);
              fill_hist_2D_parton("tt_y_tt_pt_2D", ttbar_parton.absrap(), ttbar_parton.pT()/TeV);
              if (ttbar_parton.absrap() < 0.3) fill_hist_2D_parton("tt_y_1_tt_m_t1_pt_3D", ttbar_parton.mass()/TeV, t1_parton.pT()/TeV);
              else if (ttbar_parton.absrap() < 0.9) fill_hist_2D_parton("tt_y_2_tt_m_t1_pt_3D", ttbar_parton.mass()/TeV, t1_parton.pT()/TeV);
              else if (ttbar_parton.absrap() < 2.0) fill_hist_2D_parton("tt_y_3_tt_m_t1_pt_3D", ttbar_parton.mass()/TeV, t1_parton.pT()/TeV);
            }
          }
        }

        // Get small-R jets
        const FastJets& sjets_fj = apply<FastJets>(event, "sjets");
        const Jets all_sjets = sjets_fj.jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

        // Get dressed leptons
        DressedLeptons dressedElectrons = apply<LeptonFinder>(event, "elecs").dressedLeptons();
        DressedLeptons dressedMuons     = apply<LeptonFinder>(event, "muons").dressedLeptons();

        // Perform lepton isolation
        idiscardIfAnyDeltaRLess(dressedElectrons, all_sjets, 0.4);
        idiscardIfAnyDeltaRLess(dressedMuons, all_sjets, 0.4);

        // Veto on leptons
        if (!dressedElectrons.empty()) vetoEvent;
        if (!dressedMuons.empty()) vetoEvent;

        // Get large-R jets
        const FastJets& ljets_fj = apply<FastJets>(event, "ljets");
        const Jets& trimmedJets = ljets_fj.jetsByPt();

        // Check large-R jets
        Jets ljets;
        vector<bool> b_tagged;
        for (const Jet& jet : trimmedJets) {

          if (jet.pT() < 200 * GeV)  continue;
          if (jet.pT() > 3000 * GeV) continue;
          if (jet.mass() > jet.pT()) continue;
          if (jet.mass() < 50 * GeV) continue;
          if (jet.abseta() > 2.0 )   continue;

          ljets += jet;
          b_tagged += jet.bTagged(Cuts::pT > 5 * GeV);
        }

        if (ljets.size() < 2)  vetoEvent;

        // Identify top and anti top, compute some event variables
        int top1Index(-1);
        int top2Index(-1);

        double deltaMass(FLT_MAX);
        for(int i = 0; i < (int)ljets.size(); i++) {
          if (ljets[i].pT() < 500 * GeV)  continue;
          const double diff = std::abs(ljets[i].mass() - 172.5 * GeV);
          if (diff < deltaMass) {
            deltaMass = diff;
            top1Index = i;
          }
        }

        if (top1Index == -1)  vetoEvent;

        deltaMass = FLT_MAX;
        for (int i = 0; i < (int)ljets.size(); ++i) {
          if (i == top1Index || ljets[i].pT() < 350 * GeV)  continue;
          const double diff = std::abs(ljets[i].mass() - 172.5 * GeV);
          if (diff < deltaMass) {
            deltaMass = diff;
            top2Index = i;
          }
        }

        if (top2Index == -1) vetoEvent;

        if (ljets[top1Index].pT() < ljets[top2Index].pT()) std::swap(top1Index,top2Index);

        const FourMomentum ttbar = ljets[top1Index].momentum() + ljets[top2Index].momentum();
        const FourMomentum t1 = ljets[top1Index].momentum();
        const FourMomentum t2 = ljets[top2Index].momentum();

        const double chi = calcChi(t1, t2);
        const double cosThetaStar = abs(calcCosThetaStar(t1, t2));
        if (cosThetaStar == -99) {
          MSG_DEBUG("real ttbar going faster than light! This should not happen. Vetoing event.");
          vetoEvent;
        }
        const double pout = abs(calcPout(t1, t2));
        const double dPhi = deltaPhi(t1, t2);

        // b-tagging for particle done on large-R jets
        if (!(b_tagged[top1Index] && b_tagged[top2Index]))  vetoEvent;

        // Continues with signal region cuts
        if ( abs(t1.mass() - 172.5 * GeV) > 50*GeV )  vetoEvent;
        if ( abs(t2.mass() - 172.5 * GeV) > 50*GeV )  vetoEvent;

        const FourMomentum& randomTopJet = t1; // : t2;

        _h["inclusive_particle"]->fill(0);

        fill_hist("t_pt", randomTopJet.pT()/TeV);
        fill_hist("t_y",  randomTopJet.absrap());

        fill_hist("t1_pt", t1.pT()/TeV);
        fill_hist("t1_y",  t1.absrap());
        fill_hist("t2_pt", t2.pT()/TeV);
        fill_hist("t2_y",  t2.absrap());

        fill_hist("tt_m",  ttbar.mass()/TeV);
        fill_hist("tt_pt", ttbar.pT()/TeV);
        fill_hist("tt_Ht", (t1.pT() + t2.pT())/TeV);
        fill_hist("tt_y",  ttbar.absrap());

        fill_hist("tt_yboost", 0.5 * abs(t1.rapidity() + t2.rapidity()));
        fill_hist("tt_chi", chi);
        fill_hist("tt_cosThStar", cosThetaStar);
        fill_hist("tt_pout", pout/TeV);
        fill_hist("tt_dPhi", dPhi);

        fill_hist_2D("t1_pt_t2_pt_2D", t1.pT()/TeV, t2.pT()/TeV);
        fill_hist_2D("t1_y_t2_y_2D", t1.absrap(), t2.absrap());
        fill_hist_2D("t1_y_t1_pt_2D", t1.absrap(), t1.pT()/TeV);
        fill_hist_2D("t2_y_t2_pt_2D", t2.absrap(), t2.pT()/TeV);
        fill_hist_2D("t1_pt_tt_pt_2D", t1.pT()/TeV, ttbar.pT()/TeV);
        fill_hist_2D("t1_pt_tt_m_2D", t1.pT()/TeV, ttbar.mass()/TeV);
        fill_hist_2D("tt_y_t1_pt_2D", ttbar.absrap(), t1.pT()/TeV);
        fill_hist_2D("tt_y_t1_y_2D", ttbar.absrap(), t1.absrap());
        fill_hist_2D("t1_y_tt_m_2D", t1.absrap(), ttbar.mass()/TeV);
        fill_hist_2D("tt_y_tt_m_2D", ttbar.absrap(), ttbar.mass()/TeV);
        fill_hist_2D("tt_pt_tt_m_2D", ttbar.pT()/TeV, ttbar.mass()/TeV);
        fill_hist_2D("tt_y_tt_pt_2D", ttbar.absrap(), ttbar.pT()/TeV);
        if (ttbar.absrap() < 0.3) fill_hist_2D("tt_y_1_tt_m_t1_pt_3D", ttbar.mass()/TeV, t1.pT()/TeV);
        else if (ttbar.absrap() < 0.9) fill_hist_2D("tt_y_2_tt_m_t1_pt_3D", ttbar.mass()/TeV, t1.pT()/TeV);
        else if (ttbar.absrap() < 2.0) fill_hist_2D("tt_y_3_tt_m_t1_pt_3D", ttbar.mass()/TeV, t1.pT()/TeV);

      }



      void finalize() {
        // Normalize histograms to cross-section in femtobarns (for consistency with HEPData)
        const double sf = crossSection()/picobarn * 1000 / sumOfWeights();
        for (auto& h_it : _h) {
          scale(h_it.second, sf);
          // Parton-level distributions corrected for all-hadronic branching fraction
          if (h_it.first.find("_parton") != string::npos) scale(h_it.second, 2.1882987);
          // Normalized distributions
          if (h_it.first.find("_norm") != string::npos)  normalize(h_it.second, 1.0, false);
        }
        // Multi-dimensional cross-sections
        double norm_3D = 0, norm_3D_parton = 0;
        for (auto& h_it : _h_multi) {
          if (h_it.first.find("_parton") != string::npos)  scale(h_it.second, 2.1882987);
          if (h_it.first.find("_norm") != string::npos) {
            scale(h_it.second, sf);
            if (h_it.first.find("_3D") != string::npos) {
              if (h_it.first.find("_parton") != string::npos) norm_3D_parton += h_it.second->integral(false);
              else norm_3D += h_it.second->integral(false);
              if (h_it.first.find("tt_y_1") != string::npos) { scale(h_it.second, 1. / 0.3); }
              if (h_it.first.find("tt_y_2") != string::npos) { scale(h_it.second, 1. / 0.6); }
              if (h_it.first.find("tt_y_3") != string::npos) { scale(h_it.second, 1. / 1.1); }
            }
            else {
              const double norm_2D = h_it.second->integral(false);
              scale(h_it.second, safediv(1.0, norm_2D));
              divByGroupWidth(h_it.second);
            }
          }
          else {
            if (h_it.first.find("_3D") != string::npos) {
              if (h_it.first.find("tt_y_1") != string::npos) { scale(h_it.second, 1. / 0.3); }
              if (h_it.first.find("tt_y_2") != string::npos) { scale(h_it.second, 1. / 0.6); }
              if (h_it.first.find("tt_y_3") != string::npos) { scale(h_it.second, 1. / 1.1); }
              scale(h_it.second, sf);
              divByGroupWidth(h_it.second);
            }
            else {
              scale(h_it.second, sf);
              divByGroupWidth(h_it.second);
            }
          }
        }
        scale(_h_multi["tt_y_1_tt_m_t1_pt_3D_norm"], safediv(1, norm_3D));
        divByGroupWidth(_h_multi["tt_y_1_tt_m_t1_pt_3D_norm"]);
        scale(_h_multi["tt_y_2_tt_m_t1_pt_3D_norm"], safediv(1, norm_3D));
        divByGroupWidth(_h_multi["tt_y_2_tt_m_t1_pt_3D_norm"]);
        scale(_h_multi["tt_y_3_tt_m_t1_pt_3D_norm"], safediv(1, norm_3D));
        divByGroupWidth(_h_multi["tt_y_3_tt_m_t1_pt_3D_norm"]);
        if (_mode) {
          scale(_h_multi["tt_y_1_tt_m_t1_pt_3D_parton_norm"], safediv(1, norm_3D_parton));
          divByGroupWidth(_h_multi["tt_y_1_tt_m_t1_pt_3D_parton_norm"]);
          scale(_h_multi["tt_y_2_tt_m_t1_pt_3D_parton_norm"], safediv(1, norm_3D_parton));
          divByGroupWidth(_h_multi["tt_y_2_tt_m_t1_pt_3D_parton_norm"]);
          scale(_h_multi["tt_y_3_tt_m_t1_pt_3D_parton_norm"], safediv(1, norm_3D_parton));
          divByGroupWidth(_h_multi["tt_y_3_tt_m_t1_pt_3D_parton_norm"]);
        }
      }


      double calcChi(const FourMomentum& t1, const FourMomentum& t2) {
        double ystar = 0.5 * (t1.rapidity()-t2.rapidity());
        double chi = exp( 2 * abs(ystar));
        return chi;
      }

      double calcCosThetaStar(const FourMomentum& t1, const FourMomentum& t2) {
        FourMomentum ttbar = t1 + t2;
        LorentzTransform centreOfMassTrans;
        ttbar.setX(0);
        ttbar.setY(0);
        if (ttbar.betaVec().mod2() > 1) return -99;
        centreOfMassTrans.setBetaVec( -ttbar.betaVec() );
        FourMomentum t1_star = centreOfMassTrans.transform(t1);
        double cosThetaStar;
        if (t1_star.p3().mod2() >= 0){
          cosThetaStar = t1_star.pz()/t1_star.p3().mod();
        }
        else {
          return -99;
        }
        return cosThetaStar;
      }

      double calcPout(const FourMomentum& t1, const FourMomentum& t2) {
        Vector3 t1V = t1.p3();
        Vector3 t2V = t2.p3();
        Vector3 zUnit = Vector3(0., 0., 1.);
        Vector3 vPerp = zUnit.cross(t1V);

        double pout = vPerp.dot(t2V)/vPerp.mod();
        return pout;
      }


    private:

      size_t _mode;
      map<string, Histo1DPtr> _h;
      map<string, Histo1DGroupPtr> _h_multi;

      //some functions for booking, filling and scaling the histograms
      void fill_hist(const std::string name, double value) {
        _h[name]->fill(value);
        _h[name + "_norm"]->fill(value);
      }

      void fill_hist_parton(const std::string name, double value) {
        _h[name + "_parton"]->fill(value);
        _h[name + "_parton_norm"]->fill(value);
      }

      void fill_hist_2D(const std::string name, double value_external, double value_internal) {
       _h_multi[name]->fill(value_external, value_internal);
       _h_multi[name + "_norm"]->fill(value_external, value_internal);
      }

      void fill_hist_2D_parton(const std::string name, double value_external, double value_internal) {
       _h_multi[name + "_parton"]->fill(value_external, value_internal);
       _h_multi[name + "_parton_norm"]->fill(value_external, value_internal);
      }

      void book_hist(const std::string name, unsigned int index) {
        book(_h[name], index, 1, 1);
        book(_h[name + "_norm"], index + 72, 1, 1);
        if (_mode) {
          book(_h[name + "_parton"], index + 145, 1, 1);
          book(_h[name + "_parton_norm"], index + 217, 1, 1);
        }
      }

      void book_hist_2D(const std::string name, const vector<double>& external_bins, unsigned int index) {
        book(_h_multi[name], external_bins);
        book(_h_multi[name+"_norm"], external_bins);
        for (size_t i=0; i < _h_multi[name]->numBins(); ++i) {
          book(_h_multi[name]->bin(i+1), index+i, 1, 1);
          book(_h_multi[name+"_norm"]->bin(i+1), index+72+i, 1, 1);
          if (_mode != 0) {
            book(_h_multi[name+"_parton"]->bin(i+1), index+145+i, 1, 1);
            book(_h_multi[name+"_parton_norm"]->bin(i+1), index+217+i, 1, 1);
          }
        }
      }

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2022_I2077575);
}
