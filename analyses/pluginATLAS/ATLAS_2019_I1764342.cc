// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"

namespace Rivet {


  /// @brief Z(ll)y cross-section at 13 TeV
  class ATLAS_2019_I1764342 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2019_I1764342);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // get option
      _mode = 0;
      if ( getOption("LMODE") == "EL" )  _mode = 1;
      if ( getOption("LMODE") == "MU" )  _mode = 2;

      // Prompt photons
      const PromptFinalState photon_fs(Cuts::abspid == PID::PHOTON && Cuts::pT > 30*GeV && Cuts::abseta < 2.37);
      declare(photon_fs, "Photons");

      // Prompt leptons
      const PromptFinalState bareelectron_fs = Cuts::abspid == PID::ELECTRON;
      const PromptFinalState baremuon_fs = Cuts::abspid == PID::MUON;

      // Dressed leptons
      const FinalState allphoton_fs(Cuts::abspid == PID::PHOTON); // photons used for lepton dressing
      const Cut leptoncut = Cuts::pT > 25*GeV && Cuts::abseta < 2.47;
      const LeptonFinder dressedelectron_fs(bareelectron_fs, allphoton_fs, 0.1, leptoncut);
      const LeptonFinder dressedmuon_fs(baremuon_fs, allphoton_fs, 0.1, leptoncut);

      declare(dressedelectron_fs, "Electrons");
      declare(dressedmuon_fs, "Muons");

      // FS excluding the leading photon
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(photon_fs);
      vfs.addVetoOnThisFinalState(dressedmuon_fs);
      vfs.addVetoOnThisFinalState(InvisibleFinalState());
      declare(vfs, "isolatedFS");

      // Histograms
      book(_h["EgammaT"],     2, 1, 1); // dSigma / dE^gamma_T
      book(_h["etagamma"],    3, 1, 1);
      book(_h["mZgamma"],     4, 1, 1); // dSigma / dm^{Zgamma}
      book(_h["EZgammaT"],    5, 1, 1);
      book(_h["dPhiZgamma"],  6, 1, 1);
      book(_h["ETbyMZgamma"], 7, 1, 1);
    }


   /// Perform the per-event analysis
   void analyze(const Event& event) {
     // Get objects
     Particles electrons = apply<LeptonFinder>(event, "Electrons").particlesByPt();
     Particles muons = apply<LeptonFinder>(event, "Muons").particlesByPt();
     const Particles& photons = apply<PromptFinalState>(event, "Photons").particlesByPt();

     if (photons.empty())  vetoEvent;
     if (electrons.size() < 2 && muons.size() < 2)  vetoEvent;

     if (_mode == 1 && muons.size())      vetoEvent;
     if (_mode == 2 && electrons.size())  vetoEvent;

     Particles lep;
     // Sort the dressed leptons by pt
     if (electrons.size() >= 2) {
       lep.push_back(electrons[0]);
       lep.push_back(electrons[1]);
     } else {
       lep.push_back(muons[0]);
       lep.push_back(muons[1]);
     }
     if (lep[0].Et() < 30*GeV)  vetoEvent;
     double mll = (lep[0].momentum() + lep[1].momentum()).mass();
     if (mll < 40*GeV) vetoEvent;

     Particles selectedPh;
     Particles fs = apply<VetoedFinalState>(event, "isolatedFS").particles();
     for (const Particle& ph : photons) {
       // check photon isolation
       double coneEnergy(0.0);
       for (const Particle& p : fs) {
         if ( deltaR(ph, p) < 0.2 )  coneEnergy += p.Et();
       }
       if (coneEnergy / ph.Et() > 0.07 )  continue;
       if (deltaR(ph, lep[0]) < 0.4) continue;
       if (deltaR(ph, lep[1]) < 0.4) continue;
       selectedPh.push_back(ph);
     }

     if (selectedPh.size()<1) vetoEvent;
     double mlly = (lep[0].momentum() + lep[1].momentum() + selectedPh[0].momentum()).mass();
     if (mll + mlly <= 182*GeV) vetoEvent;

     double ptlly = (lep[0].momentum() + lep[1].momentum() + selectedPh[0].momentum()).pT();
     double dphilly = deltaPhi((lep[0].momentum() + lep[1].momentum()).phi(), selectedPh[0].momentum().phi());

     // Fill plots
     _h["EgammaT"]->fill(selectedPh[0].pT()/GeV);
     _h["etagamma"]->fill(selectedPh[0].abseta());
     _h["mZgamma"]->fill(mlly/GeV);
     _h["EZgammaT"]->fill(ptlly/GeV);
     _h["dPhiZgamma"]->fill(dphilly/pi);
     _h["ETbyMZgamma"]->fill(ptlly/mlly);
   }


   /// Normalise histograms etc., after the run
   void finalize() {
      double sf = crossSection()/femtobarn/sumOfWeights();
      if (_mode == 0)  sf *= 0.5;
      scale(_h, sf);
      scale(_h["dPhiZgamma"], 1.0/pi);
   }

   /// @}


  private:

    /// Mode flag
    size_t _mode;

    /// Histograms
    map<string,Histo1DPtr> _h;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2019_I1764342);

}
