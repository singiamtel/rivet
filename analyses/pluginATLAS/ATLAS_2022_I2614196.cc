// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"

namespace Rivet {


  /// @brief Zy+jets at 13 TeV
  class ATLAS_2022_I2614196 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2022_I2614196);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Prompt photons
      const PromptFinalState photon_fs(Cuts::abspid == PID::PHOTON && Cuts::pT > 30*GeV && Cuts::abseta < 2.37);
      declare(photon_fs, "Photons");

      // Prompt leptons
      const PromptFinalState bareelectron_fs(Cuts::abspid == PID::ELECTRON);
      const PromptFinalState baremuon_fs(Cuts::abspid == PID::MUON);

      // Dressed leptons
      const FinalState allphoton_fs(Cuts::abspid == PID::PHOTON);
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

      VetoedFinalState hadrons(FinalState(Cuts::abseta < 4.4));
      hadrons.addVetoOnThisFinalState(dressedelectron_fs);
      hadrons.addVetoOnThisFinalState(dressedmuon_fs);
      declare(hadrons, "hadrons");

      //Jets
      FastJets jets(hadrons, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::DECAY);
      declare(jets, "jets");

      // Histograms
      book(_h["pTll"], 1, 1, 1);
      book(_h["DiffpTll_pTy"], 2, 1, 1);
      book(_h["SumpTll_pTy"], 3, 1, 1);
      book(_h["DeltaRll"], 4, 1, 1);
      book(_h["NJetsMix"], 5, 1, 1);
      book(_h["pTJet1"], 6, 1, 1);
      book(_h["pTJet2"], 7, 1, 1);
      book(_h["RatiopTJet12"], 8, 1, 1);
      book(_h["mjj"],  9, 1, 1);
      book(_h["mllyj"], 10, 1, 1);
      book(_h["HT"], 11, 1, 1);
      book(_h["pTysqrtHT"], 12, 1, 1);
      book(_h["DeltaPhiJetY"], 13, 1, 1);
      book(_h["pTllyj"], 14, 1, 1);
      book(_h["phi_CS_Slice1"], 51, 1, 1);
      book(_h["phi_CS_Slice2"], 51, 1, 2);
      book(_h["phi_CS_Slice3"], 51, 1, 3);
      book(_h["phi_CS_Slice4"], 51, 1, 4);
      book(_h["phi_CS_Slice5"], 51, 1, 5);
      book(_h["ctheta_CS_Slice1"], 52, 1, 1);
      book(_h["ctheta_CS_Slice2"], 52, 1, 2);
      book(_h["ctheta_CS_Slice3"], 52, 1, 3);
      book(_h["ctheta_CS_Slice4"], 52, 1, 4);
      book(_h["ctheta_CS_Slice5"], 52, 1, 5);
      book(_h["RatiopTlly_mlly_Slice1"], 17, 1, 1);
      book(_h["RatiopTlly_mlly_Slice2"], 18, 1, 1);
      book(_h["RatiopTlly_mlly_Slice3"], 19, 1, 1);
      book(_h["DiffpTll_pTy_Slice1"], 20, 1, 1);
      book(_h["DiffpTll_pTy_Slice2"], 21, 1, 1);
      book(_h["DiffpTll_pTy_Slice3"], 22, 1, 1);
      book(_h["pTllyj_Slice1"], 23, 1, 1);
      book(_h["pTllyj_Slice2"], 24, 1, 1);
      book(_h["pTllyj_Slice3"], 25, 1, 1);

      _ptllAxis = YODA::Axis<double>({0, 35, 60, 90, 135, 2500});
    }


   /// Perform the per-event analysis
   void analyze(const Event& event) {

     // Get objects
     Particles electrons = apply<LeptonFinder>(event, "Electrons").particlesByPt();
     Particles muons = apply<LeptonFinder>(event, "Muons").particlesByPt();
     Particles photons = apply<PromptFinalState>(event, "Photons").particlesByPt();

     if (photons.empty())  vetoEvent;
     if (electrons.size()<2 && muons.size()<2)  vetoEvent;

     Particles lep;
     // Sort the dressed leptons by pt
     if (electrons.size() >= 2) {
       lep.push_back(electrons[0]);
       lep.push_back(electrons[1]);
     }
     else {
       lep.push_back(muons[0]);
       lep.push_back(muons[1]);
     }

     if (lep[0].momentum().pT() < 30*GeV)  vetoEvent;

     const double mll = (lep[0].momentum() + lep[1].momentum()).mass();
     if (mll < 40*GeV)  vetoEvent;

     Particles selectedPh;
     Particles fs = apply<VetoedFinalState>(event, "isolatedFS").particles();
     for (const Particle& ph : photons) {
       // check photon isolation
       double coneEnergy(0.0);
       for (const Particle& p : fs) {
         if ( deltaR(ph, p) < 0.2 )  coneEnergy += p.Et();
       }
       if (coneEnergy / ph.momentum().pT() > 0.07 )  continue;
       if (any(electrons, deltaRLess(ph, 0.4)))  continue;
       if (any(muons, deltaRLess(ph, 0.4)))  continue;
       selectedPh.push_back(ph);
     }

     if(selectedPh.size()<1) vetoEvent;

     const double mlly   = (lep[0].momentum() + lep[1].momentum() + selectedPh[0].momentum()).mass();
     if (mll + mlly <= 182*GeV)  vetoEvent;


     // Get jets
     // abseta
     const Cut jetscut = (Cuts::pT > 30*GeV && Cuts::absrap < 2.5) || (Cuts::pT > 50*GeV && Cuts::absrapIn(2.5, 4.5));
     const Cut jetscut25 = (Cuts::pT > 25*GeV && Cuts::abseta < 2.5) || (Cuts::pT > 50*GeV && Cuts::absetaIn(2.5, 4.4));
     Jets jetsMix = apply<FastJets>(event, "jets").jetsByPt(jetscut);
     Jets jetsMix25 = apply<FastJets>(event, "jets").jetsByPt(jetscut25);
     iselectIfAnyDeltaRLess(jetsMix, photons, 0.4);
     iselectIfAnyDeltaRLess(jetsMix, electrons, 0.2);
     iselectIfAnyDeltaRLess(jetsMix25, photons, 0.4);
     iselectIfAnyDeltaRLess(jetsMix25, electrons, 0.2);

     // Jet multiplicity
     size_t njets = jetsMix.size();
     _h["NJetsMix"]->fill(njets);

     double pTJet1 = njets>0 ? jetsMix[0].pT(): -1.;
     double pTJet2 = njets>1 ? jetsMix[1].pT(): -1.;
     double mjj    = njets>1 ? (jetsMix[0].mom() + jetsMix[1].mom()).mass() : -1.;
     double mllyj  = njets>0 ? (lep[0].mom() + lep[1].mom() + selectedPh[0].mom() + jetsMix[0].mom()).mass() : -1.;

     double HT = lep[0].pT() + lep[1].pT() + selectedPh[0].pT();
     HT = sum(jetsMix25, Kin::pT, HT);

     double pTy    = selectedPh[0].pT();
     double DeltaPhijy = njets>0 ? fabs(deltaPhi(jetsMix[0], selectedPh[0])) : -1.;
     double DeltaRll   = fabs(deltaR(lep[0], lep[1]));
     double pTll   = (lep[0].momentum() + lep[1].momentum()).pT();

     double pTlly = (lep[0].mom() + lep[1].mom() + selectedPh[0].mom()).pT();
     double Ratio_pTlly_mlly = pTlly/mlly;
     double Diff_pTll_pTy    = pTll-pTy;
     double Sum_pTll_pTy     = pTll+pTy;
     double pTllyj = njets>0 ? (lep[0].mom() + lep[1].mom() + selectedPh[0].mom() + jetsMix[0].mom()).pT() : -1.;

     // Fill plots
     _h["pTJet1"]->fill(pTJet1/GeV);
     _h["pTJet2"]->fill(pTJet2/GeV);
     _h["RatiopTJet12"]->fill(pTJet2/pTJet1);
     _h["mllyj"]->fill(mllyj/GeV);
     _h["mjj"]->fill(mjj/GeV);
     _h["HT"]->fill(HT/GeV);
     _h["pTysqrtHT"]->fill((pTy/GeV)/sqrt(HT/GeV));
     _h["DeltaPhiJetY"]->fill(DeltaPhijy);
     _h["DeltaRll"]->fill(DeltaRll);
     _h["pTll"]->fill(pTll/GeV);

     if      (mlly > 125*GeV && mlly < 200*GeV)  _h["RatiopTlly_mlly_Slice1"]->fill(Ratio_pTlly_mlly);
     else if (mlly > 200*GeV && mlly < 300*GeV)  _h["RatiopTlly_mlly_Slice2"]->fill(Ratio_pTlly_mlly);
     else if (mlly > 300*GeV)                    _h["RatiopTlly_mlly_Slice3"]->fill(Ratio_pTlly_mlly);

     _h["SumpTll_pTy"]->fill(Sum_pTll_pTy/GeV);
     _h["DiffpTll_pTy"]->fill(Diff_pTll_pTy/GeV);
     if      (Sum_pTll_pTy < 200*GeV)                            _h["DiffpTll_pTy_Slice1"]->fill(Diff_pTll_pTy/GeV);
     else if (Sum_pTll_pTy > 200*GeV && Sum_pTll_pTy < 300*GeV)  _h["DiffpTll_pTy_Slice2"]->fill(Diff_pTll_pTy/GeV);
     else if (Sum_pTll_pTy > 300*GeV)                            _h["DiffpTll_pTy_Slice3"]->fill(Diff_pTll_pTy/GeV);

     _h["pTllyj"]->fill(pTllyj/GeV);
     if      (pTlly < 50*GeV)                    _h["pTllyj_Slice1"]->fill(pTllyj/GeV);
     else if (pTlly > 50*GeV && pTlly < 75*GeV)  _h["pTllyj_Slice2"]->fill(pTllyj/GeV);
     else if (pTlly > 75*GeV)                    _h["pTllyj_Slice3"]->fill(pTllyj/GeV);

     // polarization variables
     FourMomentum lep1, lep2;
     if (lep[0].pid() > 0 && lep[1].pid() < 0) {
       lep1 = lep[0].mom();
       lep2 = lep[1].mom();
     }
     else if (lep[0].pid() < 0 && lep[1].pid() > 0) {
       lep1 = lep[1].mom();
       lep2 = lep[0].mom();
     }
     else {
       MSG_DEBUG("Same sign lepton!");
       vetoEvent;
     }

     FourMomentum Dilepton = (lep1+lep2);
     FourMomentum lep1_CS;
     LorentzTransform boost;
     boost.setBetaVec(-Dilepton.betaVec());
     lep1_CS = boost.transform(lep1);
     const double ctheta_CS = cos(lep1_CS.p3().theta());
     const double phi_CS    = lep1_CS.p3().azimuthalAngle();

     size_t idx = _ptllAxis.index(pTll);
     if (mll > 80*GeV && mll<100*GeV && idx && idx <= _ptllAxis.numBins()) {
       const string s = "_Slice" + to_string(idx);
       _h["ctheta_CS" + s]->fill(ctheta_CS);
       _h["phi_CS" + s]->fill(phi_CS);
     }
   }


   /// Normalise histograms etc., after the run
   void finalize() {
     const double sf = crossSection()/femtobarn/sumOfWeights();
     scale(_h, sf);
   }

    /// @}


  private:

    /// Histograms
    map<string, Histo1DPtr> _h;

    YODA::Axis<double> _ptllAxis;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2022_I2614196);

}
