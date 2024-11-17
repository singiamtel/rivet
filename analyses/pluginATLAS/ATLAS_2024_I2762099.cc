// -*- C++ -*-

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief WZ polarisation at 13 TeV
  class ATLAS_2024_I2762099 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2024_I2762099);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Get photons to dress leptons
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);

      // Electrons and muons in Fiducial PS
      PromptFinalState leptons(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
      LeptonFinder dressedleptons(leptons, photons, 0.1, Cuts::open());
      declare(dressedleptons, "DressedLeptons");

      // Prompt neutrinos (yikes!)
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      declare(neutrinos, "Neutrinos");
      MSG_WARNING("\033[91;1mLIMITED VALIDITY - check info file for details!\033[m");

      // Muons
      PromptFinalState bare_mu(Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT);
      LeptonFinder all_dressed_mu(bare_mu, photons, 0.1, Cuts::abseta < 2.5);

      // Electrons
      PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);
      LeptonFinder all_dressed_el(bare_el, photons, 0.1, Cuts::abseta < 2.5);

      //Jet forming
      VetoedFinalState vfs(FinalState(Cuts::abseta < 5.0));
      vfs.addVetoOnThisFinalState(all_dressed_el);
      vfs.addVetoOnThisFinalState(all_dressed_mu);

      FastJets jets(vfs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::DECAY);
      declare(jets, "Jets");

      // CUT 2 PtWZ < 20 GeV
      book(_h["PtWZ20CR_DY_WZ"],  25,1,1);
      book(_h["PtWZ20CR_DY_lWZ"], 22,1,1);

      // CUT 2 PtWZ < 40 GeV
      book(_h["PtWZ40CR_DY_WZ"],  26,1,1);
      book(_h["PtWZ40CR_DY_lWZ"], 23,1,1);

      // CUT 2 PtWZ < 70 GeV
      book(_h["PtWZ70CR_DY_WZ"],  27,1,1);
      book(_h["PtWZ70CR_DY_lWZ"], 24,1,1);
    }

    void analyze(const Event& event) {

      // Selection cuts
      const vector<DressedLepton>& dressedleptons = apply<LeptonFinder>(event, "DressedLeptons").dressedLeptons();
      const Particles& neutrinos = apply<PromptFinalState>(event, "Neutrinos").particlesByPt();
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt();

      int i, j, k;
      double MassZ01 = 0., MassZ02 = 0., MassZ12 = 0.;
      double MassW0 = 0., MassW1 = 0., MassW2 = 0.;
      double WeightZ1, WeightZ2, WeightZ3;
      double WeightW1, WeightW2, WeightW3;
      double M1, M2, M3;
      double WeightTotal1, WeightTotal2, WeightTotal3;

      //---Fiducial PS: assign leptons to W and Z bosons using Resonant shape algorithm
      if (dressedleptons.size() < 3 || neutrinos.size() < 1)  vetoEvent;

      int icomb = 0;
      // try Z pair of leptons 01
      if ( (dressedleptons[0].pid() == -(dressedleptons[1].pid())) &&
           (dressedleptons[2].pid()*neutrinos[0].pid()< 0) &&
           (dressedleptons[2].abspid()==neutrinos[0].abspid()-1) ) {
        MassZ01 = (dressedleptons[0].mom() + dressedleptons[1].mom()).mass();
        MassW2 = (dressedleptons[2].mom() + neutrinos[0].mom()).mass();
        icomb = 1;
      }
      // try Z pair of leptons 02
      if ( (dressedleptons[0].pid() == -(dressedleptons[2].pid())) &&
           (dressedleptons[1].pid()*neutrinos[0].pid()< 0) &&
           (dressedleptons[1].abspid()==neutrinos[0].abspid()-1) ) {
        MassZ02 = (dressedleptons[0].mom() + dressedleptons[2].mom()).mass();
        MassW1 = (dressedleptons[1].mom() + neutrinos[0].mom()).mass();
        icomb = 2;
      }
      // try Z pair of leptons 12
      if ( (dressedleptons[1].pid() == -(dressedleptons[2].pid())) &&
         (dressedleptons[0].pid()*neutrinos[0].pid()< 0) &&
         (dressedleptons[0].abspid()==neutrinos[0].abspid()-1) ) {
        MassZ12 = (dressedleptons[1].mom() + dressedleptons[2].mom()).mass();
        MassW0 = (dressedleptons[0].mom() + neutrinos[0].mom()).mass();
        icomb = 3;
      }

      if (icomb<=0) vetoEvent;

      WeightZ1 = 1/(sqr(MassZ01*MassZ01 - MZ_PDG*MZ_PDG) + sqr(MZ_PDG*GammaZ_PDG));
      WeightW1 = 1/(sqr(MassW2*MassW2 - MW_PDG*MW_PDG) + sqr(MW_PDG*GammaW_PDG));
      WeightTotal1 = WeightZ1*WeightW1;
      M1 = -1*WeightTotal1;

      WeightZ2 = 1/(sqr(MassZ02*MassZ02- MZ_PDG*MZ_PDG) + sqr(MZ_PDG*GammaZ_PDG));
      WeightW2 = 1/(sqr(MassW1*MassW1- MW_PDG*MW_PDG) + sqr(MW_PDG*GammaW_PDG));
      WeightTotal2 = WeightZ2*WeightW2;
      M2 = -1*WeightTotal2;

      WeightZ3 = 1/(sqr(MassZ12*MassZ12 - MZ_PDG*MZ_PDG) + sqr(MZ_PDG*GammaZ_PDG));
      WeightW3 = 1/(sqr(MassW0*MassW0 - MW_PDG*MW_PDG) + sqr(MW_PDG*GammaW_PDG));
      WeightTotal3 = WeightZ3*WeightW3;
      M3 = -1*WeightTotal3;

      if( (M1 < M2 && M1 < M3) || (MassZ01 != 0 && MassW2 != 0 && MassZ02 == 0 && MassZ12 == 0) ) {
        i = 0; j = 1; k = 2;
      }
      if((M2 < M1 && M2 < M3) || (MassZ02 != 0 && MassW1 != 0 && MassZ01 == 0 && MassZ12 == 0) ) {
        i = 0; j = 2; k = 1;
      }
      if((M3 < M1 && M3 < M2) || (MassZ12 != 0 && MassW0 != 0 && MassZ01 == 0 && MassZ02 == 0) ) {
        i = 1; j = 2; k = 0;
      }

      DressedLepton Z_Lep_1 = dressedleptons[i];
      DressedLepton Z_Lep_2 = dressedleptons[j];
      DressedLepton W_Lep = dressedleptons[k];

      FourMomentum Zlepton1 = dressedleptons[i].mom();
      FourMomentum Zlepton2 = dressedleptons[j].mom();
      FourMomentum Wlepton  = dressedleptons[k].mom();
      FourMomentum Zboson   = dressedleptons[i].mom()+dressedleptons[j].mom();
      FourMomentum Wboson   = dressedleptons[k].mom()+neutrinos[0].mom();

      const double pT_WZ = (Wboson+Zboson).pT();

      double cosLepNeut;
      double norm = Wlepton.pT() * neutrinos[0].pt();
      double mTW = 0;
      if (norm != 0){
        cosLepNeut = ( Wlepton.px()*neutrinos[0].px() + Wlepton.py()*neutrinos[0].py() )/norm ;
        if(1-cosLepNeut >= 0 ) mTW = sqrt( 2 * Wlepton.pT() * neutrinos[0].pt() * (1-cosLepNeut ) );
      }

      //========================================================================================================//
      //=========================== Using RIVET to calculate all the kinematic variables =======================//
      //========================================================================================================//

      // For kinematic variable calculations
      FourMomentum WZ, W_WZframe, Z_WZframe, Wlep_WZframe, ZlepM_WZframe, ZlepP_WZframe;
      FourMomentum Wlep_Wframe, Wlep_Wframe_WZ, ZlepM_Zframe, ZlepM_Zframe_WZ, ZlepP_Zframe, ZlepP_Zframe_WZ;

      WZ = Zboson+Wboson;

      W_WZframe = Wboson;
      Z_WZframe = Zboson;
      Wlep_WZframe = Wlepton;

      FourMomentum ZleptonM, ZleptonP;

      if (Z_Lep_1.pid()>0) {
        ZleptonM = Zlepton1;
        ZleptonP = Zlepton2;
      }
      else {
        ZleptonM = Zlepton2;
        ZleptonP = Zlepton1;
      }

      ZlepM_WZframe = ZleptonM;
      ZlepP_WZframe = ZleptonP;

      LorentzTransform LT_WZ = LorentzTransform();
      LT_WZ = LT_WZ.mkObjTransformFromBeta(-(WZ.betaVec()));

      W_WZframe = LT_WZ.transform(W_WZframe);
      Z_WZframe = LT_WZ.transform(Z_WZframe);
      Wlep_WZframe = LT_WZ.transform(Wlep_WZframe);
      ZlepM_WZframe = LT_WZ.transform(ZlepM_WZframe);
      ZlepP_WZframe = LT_WZ.transform(ZlepP_WZframe);

      Wlep_Wframe = Wlepton;
      Wlep_Wframe_WZ = Wlep_WZframe;
      ZlepM_Zframe = ZleptonM;
      ZlepM_Zframe_WZ = ZlepM_WZframe;
      ZlepP_Zframe = ZleptonP;
      ZlepP_Zframe_WZ = ZlepP_WZframe;

      LorentzTransform LT_W = LorentzTransform();
      LT_W = LT_W.mkObjTransformFromBeta(-(Wboson.betaVec()));

      LorentzTransform LT_Z = LorentzTransform();
      LT_Z = LT_Z.mkObjTransformFromBeta(-(Zboson.betaVec()));

      Wlep_Wframe = LT_W.transform(Wlep_Wframe);
      ZlepM_Zframe = LT_Z.transform(ZlepM_Zframe);
      ZlepP_Zframe = LT_Z.transform(ZlepP_Zframe);

      LorentzTransform LT_W_WZ = LorentzTransform();
      LT_W_WZ = LT_W_WZ.mkObjTransformFromBeta(-(W_WZframe.betaVec()));

      LorentzTransform LT_Z_WZ = LorentzTransform();
      LT_Z_WZ = LT_Z_WZ.mkObjTransformFromBeta(-(Z_WZframe.betaVec()));

      Wlep_Wframe_WZ = LT_W_WZ.transform(Wlep_Wframe_WZ);
      ZlepM_Zframe_WZ = LT_Z_WZ.transform(ZlepM_Zframe_WZ);
      ZlepP_Zframe_WZ = LT_Z_WZ.transform(ZlepP_Zframe_WZ);

      Jets jets30;
      for (const Jet& jet : jets){
        if (jet.pt() > 30.*GeV && jet.abseta() < 4.5 &&
            deltaR(Zlepton1, jet.mom()) > 0.3 &&
            deltaR(Zlepton2, jet.mom()) > 0.3 &&
            deltaR(Wlepton, jet.mom())  > 0.3) {
          jets30 += jet;
        }
      }

      const double DY_WZ = Zboson.rap() - Wboson.rap();
      const double DY_lWZ = Zboson.rap() - Wlepton.rap();

      if (Wlepton.pT() <= 20*GeV || Zlepton1.pT() <= 15*GeV || Zlepton2.pT() <= 15*GeV)     vetoEvent;
      if (Wlepton.abseta() >= 2.5 || Zlepton1.abseta() >= 2.5 || Zlepton2.abseta() >= 2.5)  vetoEvent;
      if (fabs(Zboson.mass()/GeV - MZ_PDG) >= 10.)  vetoEvent;
      if (mTW <= 30*GeV)                            vetoEvent;
      if (deltaR(Zlepton1, Zlepton2) <= 0.2)        vetoEvent;
      if (deltaR(Zlepton1, Wlepton)  <= 0.3)        vetoEvent;
      if (deltaR(Zlepton2, Wlepton)  <= 0.3)        vetoEvent;

      // Selection cuts
      if (pT_WZ >= 70)  vetoEvent;
      _h["PtWZ70CR_DY_WZ"]->fill(fabs(DY_WZ));
      _h["PtWZ70CR_DY_lWZ"]->fill(fabs(DY_lWZ));

      if (pT_WZ >= 40)  vetoEvent;
      _h["PtWZ40CR_DY_WZ"]->fill(fabs(DY_WZ));
      _h["PtWZ40CR_DY_lWZ"]->fill(fabs(DY_lWZ));

      if (pT_WZ >= 20)  vetoEvent;
      _h["PtWZ20CR_DY_WZ"]->fill(fabs(DY_WZ));
      _h["PtWZ20CR_DY_lWZ"]->fill(fabs(DY_lWZ));

    }

    void finalize() {

      normalize(_h);

    }

    //@}


  private:


    /// @name Histograms
    //@{

    map<string,Histo1DPtr> _h;

    //@}

    const double MZ_PDG = 91.1876;
    const double MW_PDG = 83.385;
    const double GammaZ_PDG = 2.4952;
    const double GammaW_PDG = 2.085;

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2024_I2762099);
}
