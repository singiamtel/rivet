// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Tools/RivetMT2.hh"

namespace Rivet {


  /// @brief WW production in pp at 13 TeV in electroweak SUSY inspired phase space
  class ATLAS_2022_I2103950 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2022_I2103950);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs(Cuts::abseta < 5);

      // Project photons for dressing
      FinalState photons(Cuts::abspid == PID::PHOTON);

      // Cuts for leptons
      Cut lepton_cuts = ((Cuts::abseta < 2.6 && Cuts::abspid == PID::MUON ) ||
                         (Cuts::abseta < 2.47 && Cuts::abspid == PID::ELECTRON ) );

      // Muons
      PromptFinalState bare_mu(Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT);
      LeptonFinder all_dressed_mu(photons, bare_mu, 0.1);

      // Electrons
      PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);
      LeptonFinder all_dressed_el(photons, bare_el, 0.1);

      //Jet forming
      VetoedFinalState vfs(FinalState(Cuts::abseta < 5));

      InvisibleFinalState prompt_invis(OnlyPrompt::YES, TauDecaysAs::PROMPT);
      vfs.addVetoOnThisFinalState(prompt_invis);
      vfs.addVetoOnThisFinalState(all_dressed_el);
      vfs.addVetoOnThisFinalState(all_dressed_mu);

      FastJets jets(vfs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::DECAY);
      declare(jets, "jets");


      // Project dressed leptons (e/mu not from tau) with pT > 25 GeV
      PromptFinalState lep_bare(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);
      declare(lep_bare,"lep_bare");
      PromptFinalState prompt_mu(Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT);
      PromptFinalState prompt_el(Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);

      LeptonFinder lep_dressed(photons, lep_bare, 0.1, lepton_cuts, PhotonOrigin::ALL);
      declare(lep_dressed,"lep_dressed");
      LeptonFinder elecs(photons, prompt_el, 0.1, lepton_cuts, PhotonOrigin::ALL);
      declare(elecs, "elecs");
      LeptonFinder muons(photons, prompt_mu, 0.1, lepton_cuts, PhotonOrigin::ALL);
      declare(muons, "muons");

      // Get MET from generic invisibles
      VetoedFinalState ivfs(fs);
      ivfs.addVetoOnThisFinalState(VisibleFinalState(fs));
      declare(ivfs, "InvisibleFS");

      // fiducial differential cross sections
      book(_h["ptlead"], 13, 1, 1);
      book(_h["mll"], 15, 1, 1);
      book(_h["ptll"], 17, 1, 1);
      book(_h["yll"], 7, 1, 1);
      book(_h["dphill"], 9, 1, 1);
      book(_h["costhetastarll"], 11, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      const FinalState& ifs = apply<VetoedFinalState>(event, "InvisibleFS");
      const FourMomentum metinvisible = sum(ifs.particles(), FourMomentum());

      // Get met and find leptons
      //const DressedLeptons &leptons       = apply<LeptonFinder>(event, "lep_dressed").dressedLeptons();
      const DressedLeptons& all_elecs = apply<LeptonFinder>(event, "elecs").dressedLeptons();
      const DressedLeptons& all_muons = apply<LeptonFinder>(event, "muons").dressedLeptons();
      //Particles bare_leps  = apply<IdentifiedFinalState>(event, "bare_leptons").particles();

      // Find jets and jets for simplified phase space (for the latter slightly different leptons are excluded from clustering)
      Jets alljets = apply<FastJets>(event, "jets").jetsByPt(Cuts::abseta < 4.5 && Cuts::pT > 20*GeV);

      DressedLeptons elecs_muonOR;
      for (const DressedLepton& el : all_elecs) {
        bool overlaps = false;
        for (const DressedLepton& mu : all_muons) {
          if (deltaR(el, mu) < 0.01) {
            overlaps = true;
            break;
          }
        }
        if (overlaps) continue;
        elecs_muonOR.push_back(el);
      }

      // jet electron overlap removal
      for (const DressedLepton& e : elecs_muonOR) {
        ifilter_discard(alljets, deltaRLess(e, 0.2, RAPIDITY));
      }

      // muon jet overlap removal
      DressedLeptons muons;
      for (const DressedLepton& mu : all_muons) {
        float dRcut=0.4;
        if ((0.04+10.0/mu.pT()) < 0.4){
          dRcut = 0.04+10.0/mu.pT();
        }
        bool overlaps = false;
        for (const Jet& jet : alljets) {
          if (deltaR(mu, jet) < dRcut) {
            overlaps = true;
            break;
          }
        }
        if (overlaps) continue;
        muons.push_back(mu);
      }
      // electron jet overlap removal
      DressedLeptons elecs;
      for (const DressedLepton& el : elecs_muonOR) {
        float dRcut=0.4;
        if ((0.04+10/el.pT()) < 0.4){
          dRcut = 0.04+10/el.pT();
        }
        bool overlaps = false;
        for (const Jet& jet : alljets) {
          if (deltaR(el, jet) < dRcut) {
            overlaps = true;
            break;
          }
        }
        if (overlaps) continue;
        elecs.push_back(el);
      }
      Jets jets20;
      for (const Jet& jet : alljets) {
      	if (jet.abseta()<2.4 ) jets20 += jet;
      }
      // Remove events that do not contain 2 good leptons (either muons or electrons)
      if ((elecs.size()+muons.size()) !=2) vetoEvent;
      // only select electron-muon events
      if (elecs.size() != 1) vetoEvent;
      DressedLeptons leptons;
      if (elecs[0].pT()>muons[0].pT()) {
         leptons.push_back(elecs[0]);
         leptons.push_back(muons[0]);
      }
      else {
         leptons.push_back(muons[0]);
         leptons.push_back(elecs[0]);
      }


      // Define observables
      const FourMomentum dilep  = leptons.size()>1 ? leptons[0].momentum() + leptons[1].momentum() : FourMomentum(0,0,0,0);
      const double ptll         = leptons.size()>1 ? dilep.pT()/GeV : -1;
      const double Mll          = leptons.size()>1 ? dilep.mass()/GeV : -1;
      const double Yll          = leptons.size()>1 ? dilep.absrap() : -5;
      const double DPhill       = leptons.size()>1 ? fabs(deltaPhi(leptons[0], leptons[1])) : -1.;
      const double costhetastar = leptons.size()>1 ? fabs(tanh((leptons[0].eta() - leptons[1].eta()) / 2)) : -0.2;

      const double mt2 = leptons.size()>1 ? mT2(leptons[0],leptons[1], metinvisible,0) : 0;

      // Event selection for proper fiducial phase space
      if ( leptons.size() != 2)  vetoEvent;
      if ( leptons[0].pT()< 25*GeV || leptons[1].pT()< 25*GeV )  vetoEvent;

      // Veto same-flavour events
      if ( leptons[0].abspid() == leptons[1].abspid())  vetoEvent;

      // Veto same-charge events
      if ( leptons[0].pid()*leptons[1].pid()>0)  vetoEvent;

      // jet veto
      if (!(jets20.size()==0))  vetoEvent;

      if (metinvisible.pT() <= 60*GeV || metinvisible.pT() > 80*GeV )  vetoEvent;

      // m_ll cut
      if (dilep.mass() <= 100*GeV)  vetoEvent;

      //mt2 cut
      if( (mt2<=60*GeV) || (mt2>80*GeV)) vetoEvent;

      // Jet veto at 35 GeV is the default
      if (!jets20.empty())  vetoEvent;

      // fill histograms

      _h["ptlead"]->fill(leptons[0].pT()/GeV);
      _h["ptll"]->fill(ptll);
      _h["mll"]->fill(Mll);
      _h["yll"]->fill(Yll);
      _h["dphill"]->fill(DPhill);
      _h["costhetastarll"]->fill(costhetastar);

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf(crossSection()/femtobarn/sumOfWeights());
      // scale to cross section
      scale(_h, sf);
    }

    //@}

  private:

    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2022_I2103950);

}
