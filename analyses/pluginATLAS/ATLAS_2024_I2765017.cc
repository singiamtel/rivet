// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"

namespace Rivet {

  /// @brief pTmiss+jets cross-sections and ratios at 13 TeV
  class ATLAS_2024_I2765017 : public Analysis {

    public:

      /// Constructor
      RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2024_I2765017);

    public:

      // initialise
      void init() {

        _extra = getOption("MODE", "FIDUCIAL") == "EXTRAPOLATED";
        _isBSM = getOption("TYPE", "BSM") != "SM";

        ecal_cuts = Cuts::abseta < 2.47 && !(Cuts::abseta > 1.37 && Cuts::abseta < 1.52);

        // prompt photons
        PromptFinalState photons(Cuts::abspid == PID::PHOTON && Cuts::abseta < 4.9);
        declare(photons, "Photons");

        // prompt (bare) leptons
        Cut lep_inc_cuts = Cuts::abseta < 4.9 && (Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
        PromptFinalState leps(lep_inc_cuts, TauDecaysAs::PROMPT);

        // dressed leptons
        Cut el_fid_cuts  = (ecal_cuts && Cuts::abspid == PID::ELECTRON);
        Cut mu_fid_cuts  = (Cuts::abseta < 2.5 && Cuts::abspid == PID::MUON);
        Cut lep_fid_cuts = (Cuts::pT > 7*GeV) && (mu_fid_cuts || el_fid_cuts);
        LeptonFinder dressed_leps(leps, photons, 0.1, lep_fid_cuts);
        declare(dressed_leps, "DressedLeptons");

        // jet collection
        VetoedFinalState jet_fs(FinalState(Cuts::abseta < 4.9));
        if (_extra)  jet_fs.addVetoOnThisFinalState(dressed_leps);
        FastJets jets(jet_fs, JetAlg::ANTIKT, 0.4,
                      _extra? JetMuons::DECAY : JetMuons::NONE,
                      _extra? JetInvisibles::DECAY : JetInvisibles::NONE);
        declare(jets, "Jets");

        // MET
        if (_extra) declare(InvisibleFinalState(OnlyPrompt::YES), "promptNeutrinos");
        else {
          FinalState met_fs(!(Cuts::abspid == PID::MUON && (Cuts::abseta > 2.5 || Cuts::pT < 7*GeV)));
          declare(MissingMomentum(met_fs), "actualMET");
        }

        // Calorimeter particles for photon isolation
        VetoedFinalState calo_fs(VisibleFinalState(Cuts::abspid != PID::MUON));
        declare(calo_fs, "CaloParticles");

        // Jets for UE subtraction with jet-area method
        if (_extra) {
          FastJets ktjets(FinalState(), JetAlg::KT, 0.5, JetMuons::NONE, JetInvisibles::NONE);
          ktjets.useJetArea(new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec(0.9)));
          declare(ktjets, "kTjets");
        }

        // Book Histograms
        cats.init("sr0l");
        cats.init("sr0l_cr1ie_0v"); cats.init("sr0l_cr1im_0v");
        cats.init("sr0l_cr2ie_0v"); cats.init("sr0l_cr2im_0v");
        if (_extra) cats.init("sr0l_cr1ig_0v");

        for (const auto& cat : cats.regions) {
          const string mod = cat.first + "_";

          // HepData-based booking
          const string phsp1 = "met_mono";
          const string phsp2 = "met_vbf";
          const string phsp3 = "mjj_vbf";
          const string phsp4 = "dphijj_vbf";
          const string phsp4dbp = "dphij_vbf"; // Silly typo :(
          const string suff = _extra? "_dbp" : "";
          // book the ones in the paper from yoda file
          if (mod=="sr0l_") {
            book(_h[mod + phsp1], mod + phsp1 + suff);
            book(_h[mod + phsp2], mod + phsp2 + suff);
            book(_h[mod + phsp3], mod + phsp3 + suff);
            book(_h[mod + phsp4], mod + phsp4 + suff);
          } else if (mod=="sr0l_cr2im_0v_") {
            book(_h[mod + phsp1], mod + phsp1 + suff);
            book(_h[mod + phsp2], mod + phsp2 + suff);
            book(_h[mod + phsp3], mod + phsp3 + suff);
            book(_h[mod + phsp4], mod + phsp4 + suff);
          } else if (mod=="sr0l_cr1im_0v_") {
            book(_h[mod + phsp1], mod + phsp1 + suff);
            book(_h[mod + phsp2], mod + phsp2 + suff);
            book(_h[mod + phsp3], mod + phsp3 + suff);
            book(_h[mod + phsp4], mod + phsp4 + suff);
          } else if (mod=="sr0l_cr2ie_0v_") {
            book(_h[mod + phsp1], mod + phsp1 + suff);
            book(_h[mod + phsp2], mod + phsp2 + suff);
            book(_h[mod + phsp3], mod + phsp3 + suff);
            book(_h[mod + phsp4], mod + (_extra? phsp4dbp : phsp4) + suff);
          } else if (mod=="sr0l_cr1ie_0v_") {
            book(_h[mod + phsp1], mod + phsp1 + suff);
            book(_h[mod + phsp2], mod + phsp2 + suff);
            book(_h[mod + phsp3], mod + phsp3 + suff);
            book(_h[mod + phsp4], mod + phsp4 + suff);
          }
          else if (_extra && mod=="sr0l_cr1ig_0v_") {
            book(_h[mod + phsp1], mod + phsp1 + suff);
            book(_h[mod + phsp2], mod + phsp2 + suff);
            book(_h[mod + phsp3], mod + phsp3 + suff);
            book(_h[mod + phsp4], mod + phsp4 + suff);
          }

        }
      } // end of initialize

      /// Perform the per-event analysis
      void analyze(const Event& event) {

        cats.reset();

        // get prompt photons
        const Particles &photons = apply<PromptFinalState>(event, "Photons").particlesByPt(Cuts::abseta < 2.47 && Cuts::pT > 7*GeV);

        // get dressed leptons
        DressedLeptons leptons = apply<LeptonFinder>(event, "DressedLeptons").dressedLeptons(cmpMomByPt);
        // additional lepton veto
        if (leptons.size() > 2)  vetoEvent;

        // veto on *prompt* hadronic tau candidates
        for (const auto& jet : apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && ecal_cuts)) {
          for (const auto& p : jet.particles()) {
            if (p.fromHadronicTau(true))  vetoEvent; // true = only consider prompt taus
          }
        }

        // get anti-kT 0.4 jets
        Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.4);

        // remove jets within dR = 0.2 cone of a dressed lepton
        idiscardIfAnyDeltaRLess(jets, leptons, 0.2);
        // remove dressed leptons within dR = 0.4 cone of a jet
        idiscardIfAnyDeltaRLess(leptons, jets, 0.4);
        // remove jets within dR = 0.2 of a prompt photon
        if (_extra && leptons.empty() && photons.size()) {
          idiscardIfAnyDeltaRLess(jets, photons, 0.2);
        }

        const size_t njets = jets.size();
        const size_t njets_gap = centralJetVeto(jets);

        // calculate MET
        Vector3 METvec;
        if (_extra) {
          METvec = sum(apply<InvisibleFinalState>(event, "promptNeutrinos").particles(), Kin::pTvec, METvec).setZ(0);
        }
        else {
          METvec = apply<MissingMomentum>(event, "actualMET").vectorMissingPt();
        }
        const double actualMET = METvec.mod()/GeV; // actual pTmiss
        Vector3 pMETvec = sum(leptons, Kin::pTvec, METvec).setZ(0); // pseudo-MET with 'invisible' leptons
        if (_extra && leptons.empty() && photons.size()) {
          pMETvec = (photons[0].pTvec() + pMETvec).setZ(0); // pseudo-MET with 'invisible' photon
        }
        const double ptmiss = pMETvec.mod()/GeV;

        // lepton-MET kinematics
        const double lep_pt = leptons.size()? leptons[0].pT()/GeV : 0.;
        const bool hasZ = hasZcandidate(leptons);
        double dphi_metl = -1., mT_l = 0.;
        if (leptons.size()) {
          dphi_metl = deltaPhi(leptons[0], METvec); // actual MET here since SR1 and SR2 have visible leps
          mT_l = sqrt( 2 * lep_pt * actualMET * (1 - cos(dphi_metl) ) );
          if (leptons.size() > 1) {
            FourMomentum Z = leptons[0].mom() + leptons[1].mom();
          }
        }

        // photon kinematics
        const double photon_pt = photons.size()? photons[0].pT()/GeV : 0.;

        // photon isolation
        bool is_iso_gamma = true;
        if (_extra) {
          const vector<double> eta_bins = {0.0, 1.5, 3.0};
          vector<double> rho(eta_bins.size()-1, 0.0);
          FastJets kTjets = apply<FastJets>(event, "kTjets");
          for (size_t ieta = 0; ieta < eta_bins.size()-1; ++ieta) {
            fastjet::Selector fjselector(fastjet::SelectorAbsRapRange(eta_bins[ieta], eta_bins[ieta+1]));
            double sigma, area;
            kTjets.clusterSeqArea()->get_median_rho_and_sigma(fjselector, true, rho[ieta], sigma, area);
          }
          const double isoRCone = 0.4;
          for (const Particle& photon : photons) {
            // Compute calo isolation via particles within a cone around the photon
            FourMomentum mom_in_EtCone;
            for (const Particle& p : apply<VetoedFinalState>(event, "CaloParticles").particles()) {
              if (deltaR(photon, p) > isoRCone)  continue; // reject if not in cone
              mom_in_EtCone += p; // sum momentum
            }
            mom_in_EtCone -= photon; // subtract core photon
            // UE subtraction energy
            const double UEpT = M_PI * sqr(isoRCone) * rho[binIndex(photon.abseta(), eta_bins)];
            // Use photon if energy in isolation cone is low enough
            if (mom_in_EtCone.Et()/GeV - UEpT > 0.044*photon.mom().pT()/GeV + 2.45) is_iso_gamma = false;
            break;
          }
        }

        // jet kinematics
        double jpt1 = 0, jeta1 = 0., jpt2 = 0.;
        double mjj = 0., drapjj = 0., dphijj = 0.;
        //double dphi_metj = 0., mT_j = 0.;
        if (njets) {
          jpt1 = jets[0].pT()/GeV;
          jeta1 = jets[0].eta();
          if (njets >= 2) {
            mjj = (jets[0].mom() + jets[1].mom()).mass()/GeV;
            jpt2 = jets[1].pT()/GeV;
            drapjj = deltaRap(jets[0], jets[1]);
            dphijj = signedDeltaPhi(jets[0], jets[1]);
          }
        }

        // jet-MET balance (check dPhi between MET and first 4 jets)
        bool fail_dphi_jet_MET = false;
        for (size_t i = 0; i < jets.size() && i < 4; ++i) {
          fail_dphi_jet_MET |= (deltaPhi(jets[i], pMETvec) < 0.4);
        }

        // start categorising
        if (leptons.empty() && ptmiss > 200. && !fail_dphi_jet_MET) { // 0-lep SRs
          cats.trigger("sr0l");
        }

        if (leptons.size() == 1 && ptmiss > 200. && !fail_dphi_jet_MET) { // CR for 0-lep SR using W
          // higher leading electron pT cut due to trigger, cut on actual MET to suppress fakes
          if (leptons[0].abspid() == PID::ELECTRON && lep_pt > 30. && actualMET > 60. && inRange(mT_l, 30., 100.)) {
            cats.trigger("sr0l_cr1ie_0v");
          }
          else if (leptons[0].abspid() == PID::MUON) {
            cats.trigger("sr0l_cr1im_0v");
          }
        }
        if (hasZ && ptmiss > 200. && lep_pt > 80. && !fail_dphi_jet_MET) { // CR for 0-lep SR using Z
          if (leptons[0].abspid() == PID::ELECTRON) {
            cats.trigger("sr0l_cr2ie_0v");
          }
          else if (leptons[0].abspid() == PID::MUON) {
            cats.trigger("sr0l_cr2im_0v");
          }
        }
        if (_extra && leptons.empty() && ptmiss > 200. && photons.size() && photon_pt > 160.) {
          if (!fail_dphi_jet_MET && is_iso_gamma) { // CR for 0-lep SR using gamma
            cats.trigger("sr0l_cr1ig_0v");
          }
        }

        // identify jet topology
        const bool pass_monojet = jpt1 > 120*GeV && fabs(jeta1) < 2.4;
        const bool pass_vbf = njets >= 2 && mjj > 200*GeV && jpt1 > 80. && jpt2 > 50. && fabs(drapjj) > 1.; // CJV applied below

        // fill histograms for all categories
        for (const auto& cat : cats.regions) {

          // check if category was triggered
          if (!cat.second)  continue;
          // construct prefix for histogram name
          const string mod = cat.first + "_";

          if (pass_monojet) {
            _h[mod + "met_mono"]->fill(ptmiss);
          } // end of pass monojet

          if (pass_vbf) {
            if (!njets_gap) { // gap-jet veto
              // This is the VBF region!
              _h[mod + "met_vbf"]->fill(ptmiss);
              _h[mod + "mjj_vbf"]->fill(mjj);
              _h[mod + "dphijj_vbf"]->fill(dphijj/pi);
            }
          } // end of pass VBF

        } // end of loop over categories

      }// end of analyze


      /// Normalise, scale and otherwise manipulate histograms here
      void finalize() {
        scale(_h, crossSection() / sumOfWeights() / femtobarn);

        if (_isBSM && !_extra) {
          for (const string& auxil : vector<string>{"cr1ie_"s, "cr1im_"s, "cr2ie_"s, "cr2im_"s}) {
            for (const string& obs : vector<string>{"met_mono"s, "met_vbf"s, "mjj_vbf"s, "dphijj_vbf"s}) {
              const string rmisslabel("rmiss_"s+auxil+obs);
              const string statslabel("stats");
              MSG_DEBUG("Constructing Rmiss for " << rmisslabel);
              if (_e.find(rmisslabel) == _e.end())  book(_e[rmisslabel], rmisslabel); // book for first weight
              // load SM prediction
              const string numlabel("sr0l_"s+obs);
              const string denlabel("sr0l_"s+auxil+"0v_"s+obs);
              const YODA::Estimate1D& numSM = refData(numlabel+"_thy_nlo"s);
              const YODA::Estimate1D& denSM = refData(denlabel+"_thy_nlo"s);
              const YODA::Estimate1D& rmiss = refData(rmisslabel+"_thy_nlo"s);
              // inject BSM component
              YODA::Estimate1D numer = numSM + _h[numlabel]->mkEstimate("", statslabel);
              YODA::Estimate1D denom = denSM + _h[denlabel]->mkEstimate("", statslabel);
              // construct Rmiss
              if (numer != denom)  numer.rebinXTo(denom.xEdges()); // harmonise binning if need be
              YODA::Estimate1D ratio = numer / denom; // assumes uncorrelated error breakdowns
              // copy over old Rmiss error breakdown
              const string path = _e[rmisslabel]->path();
              *_e[rmisslabel] = rmiss; _e[rmisslabel]->setPath(path);
              // update central value and stats error component
              for (auto& b : _e[rmisslabel]->bins()) {
                const auto& injb = ratio.bin(b.index());
                b.set(injb.val(), injb.err(statslabel), statslabel);
              }
            }
          }
        }
      }


      // check if jet is between tagging jets
      bool isBetween(const Jet& probe, const Jet& boundary1, const Jet& boundary2) const {
        double y_p = probe.rapidity();
        double y_b1 = boundary1.rapidity();
        double y_b2 = boundary2.rapidity();

        double y_min = std::min(y_b1, y_b2);
        double y_max = std::max(y_b1, y_b2);

        return  (y_p > y_min && y_p < y_max);
      }


      // count number of gap jets for central jet veto
      size_t centralJetVeto(const Jets &jets) const {
        if (jets.size() < 2) return 0;
        const Jet& bj1 = jets.at(0);
        const Jet& bj2 = jets.at(1);

        size_t n_between = 0;
        // start loop at the 3rd hardest pT jet
        for (size_t i = 2; i < jets.size(); ++i) {
          const Jet& j = jets.at(i);
          if (isBetween(j, bj1, bj2))  ++n_between;
        }
        return n_between;
      }


      bool hasZcandidate(const DressedLeptons& leptons) const {
        if (leptons.size() != 2)  return false; // ask for exactly two leptons
        if (leptons[0].pid() != -leptons[1].pid())  return false; // check for SFOS
        double Zmass = (leptons[0].mom() + leptons[1].mom()).mass();
        return inRange(Zmass, 66*GeV,  116*GeV); // check dilepton mass
      }


      double signedDeltaPhi(const Jet& j1, const Jet& j2) const {
        double dphijj;
        if (j1.rap() > j2.rap()) {
          dphijj = j1.phi() - j2.phi();
        }
        else {
          dphijj = j2.phi() - j1.phi();
        }
        return mapAngleMPiToPi(dphijj);
      }


      struct Categories {
        map<string, bool> regions;

        Categories () { }

        void init(const string& name) { regions[name] = false; }

        void trigger(const string& name) {
          regions[name] = true;
        }

        void reset() {
          for (auto& reg : regions) { reg.second = false; }
        }
      };


    private:

      // Data members like post-cuts event weight counters go here
      /// @name Histograms
      //@{
      map<string, Histo1DPtr> _h;

      map<string, Estimate1DPtr> _e;

      Categories cats;

      Cut ecal_cuts;

      size_t _extra, _isBSM;
      //@}

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2024_I2765017);
}

