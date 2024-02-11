// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// @brief MC validation analysis for inclusive W events
  class MC_WINC : public Analysis {
  public:

    /// Default constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_WINC);


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {

      // Use analysis options
      _dR = (getOption("SCHEME") == "BARE") ? 0.0 : 0.2;
      _lepton = (getOption("LMODE") == "MU") ? PID::MUON : PID::ELECTRON;
      const double ETACUT = getOption<double>("ABSETALMAX", 3.5);
      const double PTCUT = getOption<double>("PTLMIN", 25.);
      const Cut cut = Cuts::abseta < ETACUT && Cuts::pT > PTCUT*GeV;

      // Define projections
      declare("MET", MissingMomentum());
      LeptonFinder lf(_dR, cut && Cuts::abspid == _lepton);
      declare(lf, "Leptons");

      double sqrts = sqrtS() > 0.0 ? sqrtS() : 14000.;
      book(_h_W_mass ,"W_mass", 50, 55.0, 105.0);
      book(_h_W_mT ,"W_mT", 40, 60.0, 100.0);
      book(_h_W_pT ,"W_pT", logspace(100, 1.0, 0.5*sqrts));
      book(_h_W_pT_peak ,"W_pT_peak", 25, 0.0, 125.0);
      book(_h_W_y ,"W_y", 40, -4.0, 4.0);
      book(_h_W_phi ,"W_phi", 25, 0.0, TWOPI);
      book(_h_Wplus_pT ,"Wplus_pT", logspace(25, 10.0, 0.5*sqrts));
      book(_h_Wminus_pT ,"Wminus_pT", logspace(25, 10.0, 0.5*sqrts));
      book(_h_lepton_pT ,"lepton_pT", logspace(100, 10.0, 0.25*sqrts));
      book(_h_lepton_eta ,"lepton_eta", 40, -4.0, 4.0);
      book(_htmp_dsigminus_deta ,"lepton_dsigminus_deta", 20, 0.0, 4.0);
      book(_htmp_dsigplus_deta  ,"lepton_dsigplus_deta", 20, 0.0, 4.0);

      book(_h_asym, "W_chargeasymm_eta");
      book(_h_asym_pT, "W_chargeasymm_pT");
    }



    /// Do the analysis
    void analyze(const Event & event) {

      // MET cut
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      if (pmiss.pT() < 25*GeV) vetoEvent;

      // Identify the closest-matching l+MET to m == mW
      const Particles& ls = apply<LeptonFinder>(event, "Leptons").particles();
      const int ifound = closestMatchIndex(ls, pmiss, Kin::mass, 80.4*GeV, 60*GeV, 100*GeV);
      if (ifound < 0) vetoEvent;
      const Particle& l = ls[ifound];

      double charge3_x_eta = 0;
      int charge3 = 0;
      FourMomentum wmom = l.mom() + pmiss;
      _h_W_mass->fill(wmom.mass()/GeV);
      _h_W_mT->fill(mT(l, pmiss)/GeV);
      _h_W_pT->fill(wmom.pT()/GeV);
      _h_W_pT_peak->fill(wmom.pT()/GeV);
      _h_W_y->fill(wmom.rapidity());
      _h_W_phi->fill(wmom.phi());
      _h_lepton_pT->fill(l.pT()/GeV);
      _h_lepton_eta->fill(l.eta());
      FourMomentum emom;
      if (l.charge3() != 0) {
        emom = l.momentum();
        charge3_x_eta = l.charge3() * emom.eta();
        charge3 = l.charge3();
      }
      assert(charge3_x_eta != 0);
      assert(charge3!=0);
      if (emom.Et() > 30/GeV) {
        if (charge3_x_eta < 0) {
          _htmp_dsigminus_deta->fill(emom.eta());
        } else {
          _htmp_dsigplus_deta->fill(emom.eta());
        }
      }
      if (charge3 < 0) {
        _h_Wminus_pT->fill(wmom.pT()/GeV);
      } else {
        _h_Wplus_pT->fill(wmom.pT()/GeV);
      }
    }


    /// Finalize
    void finalize() {
      scale(_h_W_mass, crossSection()/picobarn/sumOfWeights());
      scale(_h_W_mT, crossSection()/picobarn/sumOfWeights());
      scale(_h_W_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_W_pT_peak, crossSection()/picobarn/sumOfWeights());
      scale(_h_W_y, crossSection()/picobarn/sumOfWeights());
      scale(_h_W_phi, crossSection()/picobarn/sumOfWeights());
      scale(_h_lepton_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_lepton_eta, crossSection()/picobarn/sumOfWeights());

      // Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta) for each Et region
      divide(*_htmp_dsigplus_deta - *_htmp_dsigminus_deta,
             *_htmp_dsigplus_deta + *_htmp_dsigminus_deta,
             _h_asym);

      // // W charge asymmetry vs. pTW: dsig+/dpT / dsig-/dpT
      divide(_h_Wplus_pT, _h_Wminus_pT,
             _h_asym_pT);

      scale(_h_Wplus_pT, crossSection()/picobarn/picobarn/sumOfWeights());
      scale(_h_Wminus_pT, crossSection()/picobarn/picobarn/sumOfWeights());
    }

    /// @}


  protected:

    /// @name Parameters for specialised e/mu and dressed/bare subclassing
    /// @{
    double _dR;
    PdgId _lepton;
    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_W_mass;
    Histo1DPtr _h_W_mT;
    Histo1DPtr _h_W_pT;
    Histo1DPtr _h_W_pT_peak;
    Histo1DPtr _h_W_y;
    Histo1DPtr _h_W_phi;
    Histo1DPtr _h_Wplus_pT;
    Histo1DPtr _h_Wminus_pT;
    Histo1DPtr _h_lepton_pT;
    Histo1DPtr _h_lepton_eta;

    Histo1DPtr _htmp_dsigminus_deta;
    Histo1DPtr _htmp_dsigplus_deta;

    Estimate1DPtr _h_asym;
    Estimate1DPtr _h_asym_pT;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(MC_WINC);

}
