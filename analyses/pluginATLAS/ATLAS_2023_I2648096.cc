#include "Rivet/Analysis.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/HistoGroup.hh"

namespace Rivet {

  /// @brief lepton differential ttbar analysis at 13 TeV
  class ATLAS_2023_I2648096 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2023_I2648096);

    void init() {

      Cut eta_full = Cuts::abseta < 5.0 && Cuts::pT > 1.0*MeV;

      // Get photons to dress leptons
      FinalState photons(Cuts::abspid == PID::PHOTON);

      // Projection to find the electrons
      PromptFinalState prompt_el(Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);
      LeptonFinder elecs(prompt_el, photons, 0.1, (Cuts::abseta < 2.5 && Cuts::pT > 25*GeV));
      LeptonFinder veto_elecs(prompt_el, photons, 0.1, eta_full);
      declare(elecs, "elecs");

      // Projection to find the muons
      PromptFinalState prompt_mu(Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT);
      LeptonFinder muons(prompt_mu, photons, 0.1, (Cuts::abseta < 2.5 && Cuts::pT > 25*GeV));
      LeptonFinder veto_muons(prompt_mu, photons, 0.1, eta_full);
      declare(muons, "muons");

      const InvisibleFinalState invis(OnlyPrompt::YES, TauDecaysAs::PROMPT);
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(veto_elecs);
      vfs.addVetoOnThisFinalState(veto_muons);
      vfs.addVetoOnThisFinalState(invis);
      FastJets jets(vfs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::ALL);
      declare(jets, "jets");

      // Book histograms
      bookHistos("lep_pt",       6);
      bookHistos("lep_eta",      9);
      bookHistos("dilep_sumE",  12);
      bookHistos("dilep_mass",  15);
      bookHistos("dilep_sumpt", 18);
      bookHistos("dilep_pt",    21);
      bookHistos("dilep_dphi",  24);
      bookHistos("dilep_rap",   27);

      // unrolled 2D distributions - 2nd-dim bin edges must be specified
      const vector<double> b_mass_2D = { 0., 70., 100., 130., 200., 800. };
      const vector<double> b_pt_2D = { 0., 40., 65., 100. };
      const vector<double> b_sumE_2D = { 0., 110., 140., 200., 250., 900. };

      bookHisto2D("dilep_rap_mass",  78, b_mass_2D);
      bookHisto2D("dilep_dphi_mass", 79, b_mass_2D);
      bookHisto2D("dilep_dphi_pt",   80, b_pt_2D);
      bookHisto2D("dilep_dphi_sumE", 81, b_sumE_2D);
    }

    void analyze(const Event& event) {
      DressedLeptons elecs = apply<LeptonFinder>(event, "elecs").dressedLeptons();
      DressedLeptons muons = apply<LeptonFinder>(event, "muons").dressedLeptons();
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

      // Check overlap of jets/leptons.
      for (const Jet& jet : jets) {
        idiscard(elecs, deltaRLess(jet, 0.4));
        idiscard(muons, deltaRLess(jet, 0.4));
      }
      if (elecs.empty() || muons.empty())  vetoEvent;
      if (elecs[0].charge() == muons[0].charge())  vetoEvent;

      FourMomentum el = elecs[0].momentum();
      FourMomentum mu = muons[0].momentum();
      FourMomentum ll = elecs[0].momentum() + muons[0].momentum();

      if (max(el.pT(), mu.pT()) < 27*GeV)  vetoEvent;

      const double dphi = deltaPhi(el,mu);
      const double drap = ll.absrap();
      const double pt   = ll.pT()/GeV;
      const double sumE = min((el.E()+mu.E())/GeV, 899);
      const double mass = min(ll.mass()/GeV, 799.);

      // Fill histograms
      // include explicit overflow protection as last bins are inclusive
      fillHistos("lep_pt",      min(el.pT()/GeV, 349.));
      fillHistos("lep_pt",      min(mu.pT()/GeV, 349.));
      fillHistos("lep_eta",     el.abseta());
      fillHistos("lep_eta",     mu.abseta());
      fillHistos("dilep_pt",    min(pt, 299.));
      fillHistos("dilep_mass",  mass);
      fillHistos("dilep_rap",   drap);
      fillHistos("dilep_dphi",  dphi);
      fillHistos("dilep_sumpt", min((el.pT()+mu.pT())/GeV, 599));
      fillHistos("dilep_sumE",  sumE);

      // Fill unrolled 2D histograms vs mass
      fillHisto2D("dilep_rap_mass",  mass, drap);
      fillHisto2D("dilep_dphi_mass", mass, dphi);
      fillHisto2D("dilep_dphi_pt",   min(pt, 99.), dphi);
      fillHisto2D("dilep_dphi_sumE", sumE, dphi);
    }

    void finalize() {
      // Normalize to cross-section
      const double sf = crossSection()/femtobarn/sumOfWeights();

      // finalisation of 1D histograms
      for (auto& hist : _h) {
        const double norm = 1.0 / hist.second->integral();
        // histogram normalisation
        if (hist.first.find("norm") != string::npos)  scale(hist.second, norm);
        else  scale(hist.second, sf);
      }

      // finalisation of 2D histograms
      for (auto& hist : _h_multi) {
        if (hist.first.find("_norm") != std::string::npos) {
          const double sf = safediv(1.0, hist.second->integral(false));
          scale(hist.second, sf);
        }
        else {
          scale(hist.second, sf);
        }
        divByGroupWidth(hist.second);
      }
    }


  private:

    /// @name Histogram helper functions
    //@{
    void bookHistos(const string& name, const unsigned int id) {
      book(_h[name], id, 1, 1);
      book(_h["norm_" + name], id+36, 1, 1);
    }

    void fillHistos(const string& name, const double value) {
      _h[name]->fill(value);
      _h["norm_" + name]->fill(value);
    }

    void bookHisto2D(const string& name, const unsigned int id, const vector<double>& axis1) {
      book(_h_multi[name], axis1);
      book(_h_multi[name + "_norm"], axis1);
      for (size_t i = 1; i < axis1.size(); ++i) {
        book(_h_multi[name]->bin(i), id, 1, i);
        book(_h_multi[name + "_norm"]->bin(i), id+4, 1, i);
      }
    }



    void fillHisto2D(const string& name, const double val1, const double val2) {
      _h_multi[name]->fill(val1, val2);
      _h_multi[name+"_norm"]->fill(val1, val2);
    }


    // pointers to 1D and 2D histograms
    map<string, Histo1DPtr> _h;
    map<string, Histo1DGroupPtr> _h_multi;
    //@}
    // acceptance counter

  };

  // Declare the class as a hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2023_I2648096);
}
