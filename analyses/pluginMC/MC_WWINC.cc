// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief MC validation analysis for W^+[enu]W^-[munu] events
  class MC_WWINC : public Analysis {
  public:

    /// Default constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_WWINC);


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      declare("MET", MissingMomentum());

      // Find electrons with cuts from input options
      const double ETAECUT = getOption<double>("ABSETAEMAX", 3.5);
      const double PTECUT = getOption<double>("PTEMIN", 25.);
      const Cut cut_e = Cuts::abseta < ETAECUT && Cuts::pT > PTECUT*GeV;
      LeptonFinder ef(0.2, cut_e && Cuts::abspid == PID::ELECTRON);
      declare(ef, "Elecs");

      // Find muons with cuts from input options
      const double ETAMUCUT = getOption<double>("ABSETAMUMAX", 3.5);
      const double PTMUCUT = getOption<double>("PTMUMIN", 25.);
      const Cut cut_m = Cuts::abseta < ETAMUCUT && Cuts::pT > PTMUCUT*GeV;
      LeptonFinder mf(0.2, cut_m && Cuts::abspid == PID::MUON);
      declare(mf, "Muons");


      // properties of the pair momentum
      double sqrts = sqrtS()>0. ? sqrtS() : 14000.;
      book(_h_WW_pT ,"WW_pT", logspace(100, 1.0, max(1.1,0.5*sqrts)));
      book(_h_WW_pT_peak ,"WW_pT_peak", 25, 0.0, 25.0);
      book(_h_WW_eta ,"WW_eta", 40, -7.0, 7.0);
      book(_h_WW_phi ,"WW_phi", 25, 0.0, TWOPI);
      book(_h_WW_m ,"WW_m", logspace(100, 150.0, 180.0+0.25*sqrts));

      // correlations between the WW
      book(_h_WW_dphi ,"WW_dphi", 25, 0.0, PI);  /// @todo non-linear?
      book(_h_WW_deta ,"WW_deta", 25, -7.0, 7.0);
      book(_h_WW_dR ,"WW_dR", 25, 0.5, 7.0);
      book(_h_WW_dpT ,"WW_dpT", logspace(100, 1.0, max(1.1,0.5*sqrts)));
      book(_h_WW_costheta_planes ,"WW_costheta_planes", 25, -1.0, 1.0);

      /// @todo fuer WW: missing ET

      // properties of the W bosons
      book(_h_W_pT ,"W_pT", logspace(100, 10.0, max(11.,0.25*sqrts)));
      book(_h_W_eta ,"W_eta", 70, -7.0, 7.0);

      // properties of the leptons
      book(_h_Wl_pT ,"Wl_pT", logspace(100, 30.0, max(31., 0.1*sqrts)));
      book(_h_Wl_eta ,"Wl_eta", 40, -3.5, 3.5);

      // correlations between the opposite charge leptons
      book(_h_WeWm_dphi ,"WeWm_dphi", 25, 0.0, PI);
      book(_h_WeWm_deta ,"WeWm_deta", 25, -5.0, 5.0);
      book(_h_WeWm_dR ,"WeWm_dR", 25, 0.5, 5.0);
      book(_h_WeWm_m ,"WeWm_m", 100, 0.0, 300.0);
    }


    /// Do the analysis
    void analyze(const Event& event) {

      // MET cut
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      if (pmiss.pT() < 25*GeV) vetoEvent;

      // Identify the closest-matching l+MET to m == mW
      /// @note Dubious strategy, given there are two neutrinos...
      const Particles& es = apply<LeptonFinder>(event, "Elecs").particles();
      const int iefound = closestMatchIndex(es, pmiss, Kin::mass, 80.4*GeV, 60*GeV, 100*GeV);
      const Particles& mus = apply<LeptonFinder>(event, "Muons").particles();
      const int imfound = closestMatchIndex(mus, pmiss, Kin::mass, 80.4*GeV, 60*GeV, 100*GeV);

      // Require two valid W candidates
      if (iefound < 0 || imfound < 0) vetoEvent;

      // Get momenta
      const FourMomentum pe = es[iefound].mom();
      const FourMomentum pm = mus[imfound].mom();

      const FourMomentum pww = pe + pm + pmiss; //< don't double-count the MET
      _h_WW_pT->fill(pww.pT());
      _h_WW_pT_peak->fill(pww.pT());
      _h_WW_eta->fill(pww.eta());
      _h_WW_phi->fill(pww.phi());
      if (pww.mass2() > 0.0) _h_WW_m->fill(pww.mass());

      const FourMomentum penu = pe + pmiss;
      const FourMomentum pmnu = pm + pmiss;
      _h_WW_dphi->fill(mapAngle0ToPi(penu.phi()-pmnu.phi()));
      _h_WW_deta->fill(penu.eta() - pmnu.eta());
      _h_WW_dR->fill(deltaR(penu, pmnu));
      _h_WW_dpT->fill(fabs(penu.pT() - pmnu.pT()));

      const Vector3 crossWenu = pe.p3().cross(penu.p3());
      const Vector3 crossWmnu = pm.p3().cross(pmnu.p3());
      const double costheta = crossWenu.dot(crossWmnu)/crossWenu.mod()/crossWmnu.mod();
      _h_WW_costheta_planes->fill(costheta);

      _h_W_pT->fill(penu.pT()/GeV);
      _h_W_pT->fill(pmnu.pT()/GeV);
      _h_W_eta->fill(penu.eta());
      _h_W_eta->fill(pmnu.eta());

      _h_Wl_pT->fill(pe.pT()/GeV);
      _h_Wl_pT->fill(pm.pT()/GeV);
      _h_Wl_eta->fill(pe.eta());
      _h_Wl_eta->fill(pm.eta());

      _h_WeWm_dphi->fill(mapAngle0ToPi(pe.phi() - pm.phi()));
      _h_WeWm_deta->fill(pe.eta() - pm.eta());
      _h_WeWm_dR->fill(deltaR(pe, pm));
      const double m2 = max(FourMomentum(pe + pm).mass2(), 0.0);
      _h_WeWm_m->fill(sqrt(m2));
    }


    /// Finalize
    void finalize() {
      const double norm = crossSection()/picobarn/sumOfWeights();
      scale(_h_WW_pT, norm);
      scale(_h_WW_pT_peak, norm);
      scale(_h_WW_eta, norm);
      scale(_h_WW_phi, norm);
      scale(_h_WW_m, norm);
      scale(_h_WW_dphi, norm);
      scale(_h_WW_deta, norm);
      scale(_h_WW_dR, norm);
      scale(_h_WW_dpT, norm);
      scale(_h_WW_costheta_planes, norm);
      scale(_h_W_pT, norm);
      scale(_h_W_eta, norm);
      scale(_h_Wl_pT, norm);
      scale(_h_Wl_eta, norm);
      scale(_h_WeWm_dphi, norm);
      scale(_h_WeWm_deta, norm);
      scale(_h_WeWm_dR, norm);
      scale(_h_WeWm_m, norm);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_WW_pT;
    Histo1DPtr _h_WW_pT_peak;
    Histo1DPtr _h_WW_eta;
    Histo1DPtr _h_WW_phi;
    Histo1DPtr _h_WW_m;
    Histo1DPtr _h_WW_dphi;
    Histo1DPtr _h_WW_deta;
    Histo1DPtr _h_WW_dR;
    Histo1DPtr _h_WW_dpT;
    Histo1DPtr _h_WW_costheta_planes;
    Histo1DPtr _h_W_pT;
    Histo1DPtr _h_W_eta;
    Histo1DPtr _h_Wl_pT;
    Histo1DPtr _h_Wl_eta;
    Histo1DPtr _h_WeWm_dphi;
    Histo1DPtr _h_WeWm_deta;
    Histo1DPtr _h_WeWm_dR;
    Histo1DPtr _h_WeWm_m;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(MC_WWINC);

}
