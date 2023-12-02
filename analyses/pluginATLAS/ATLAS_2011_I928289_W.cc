// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeptonFinder.hh"

namespace Rivet {


  /// W inclusive cross-sections at 7 TeV
  class ATLAS_2011_I928289_W : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_I928289_W);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare("MET", MissingMomentum());
      LeptonFinder ef_bare(0.0, Cuts::pT > 20*GeV && Cuts::abspid == PID::ELECTRON);
      LeptonFinder ef_dressed(0.1, Cuts::pT > 20*GeV && Cuts::abspid == PID::ELECTRON);
      LeptonFinder mf_bare(0.0, Cuts::pT > 20*GeV && Cuts::abspid == PID::MUON);
      LeptonFinder mf_dressed(0.1, Cuts::pT > 20*GeV && Cuts::abspid == PID::MUON);
      declare(ef_bare, "Elecs_bare");
      declare(ef_dressed, "Elecs_dressed");
      declare(mf_bare, "Muons_bare");
      declare(mf_dressed, "Muons_dressed");

      /// Book histograms here
      book(_h_Wminus_lepton_eta_el_bare   ,3, 1, 1);
      book(_h_Wminus_lepton_eta_el_dressed,3, 1, 2);
      book(_h_Wminus_lepton_eta_mu_bare   ,3, 1, 3);
      book(_h_Wminus_lepton_eta_mu_dressed,3, 1, 4);
      book(_h_Wplus_lepton_eta_el_bare    ,5, 1, 1);
      book(_h_Wplus_lepton_eta_el_dressed ,5, 1, 2);
      book(_h_Wplus_lepton_eta_mu_bare    ,5, 1, 3);
      book(_h_Wplus_lepton_eta_mu_dressed ,5, 1, 4);
      book(_h_W_asym_eta_el_bare          ,7, 1, 1);
      book(_h_W_asym_eta_el_dressed       ,7, 1, 2);
      book(_h_W_asym_eta_mu_bare          ,7, 1, 3);
      book(_h_W_asym_eta_mu_dressed       ,7, 1, 4);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // W reco, starting with MET
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      if (pmiss.Et() < 25*GeV) vetoEvent;

      const Particles& esb = apply<LeptonFinder>(event, "Elecs_bare").particles();
      const Particles esb_mtfilt = select(esb, [&](const Particle& e){ return mT(e, pmiss) > 40*GeV; });
      const int iebfound = closestMatchIndex(esb_mtfilt, pmiss, Kin::mass, 80.4*GeV, 40*GeV);
      if (iebfound >= 0) {
        const Particle& l = esb_mtfilt[iebfound];
        (l.charge3() > 0 ? _h_Wplus_lepton_eta_el_bare : _h_Wminus_lepton_eta_el_bare)->fill(l.abseta());
      }

      const Particles& esd = apply<LeptonFinder>(event, "Elecs_dressed").particles();
      const Particles esd_mtfilt = select(esd, [&](const Particle& e){ return mT(e, pmiss) > 40*GeV; });
      const int iedfound = closestMatchIndex(esd_mtfilt, pmiss, Kin::mass, 80.4*GeV, 40*GeV);
      if (iedfound >= 0) {
        const Particle& l = esd_mtfilt[iedfound];
        (l.charge3() > 0 ? _h_Wplus_lepton_eta_el_dressed : _h_Wminus_lepton_eta_el_dressed)->fill(l.abseta());
      }

      const Particles& musb = apply<LeptonFinder>(event, "Muons_bare").particles();
      const Particles musb_mtfilt = select(musb, [&](const Particle& m){ return mT(m, pmiss) > 40*GeV; });
      const int imbfound = closestMatchIndex(musb_mtfilt, pmiss, Kin::mass, 80.4*GeV, 40*GeV);
      if (imbfound >= 0) {
        const Particle& l = musb_mtfilt[imbfound];
        (l.charge3() > 0 ? _h_Wplus_lepton_eta_mu_bare : _h_Wminus_lepton_eta_mu_bare)->fill(l.abseta());
      }

      const Particles& musd = apply<LeptonFinder>(event, "Muons_dressed").particles();
      const Particles musd_mtfilt = select(musd, [&](const Particle& m){ return mT(m, pmiss) > 40*GeV; });
      const int imdfound = closestMatchIndex(musd_mtfilt, pmiss, Kin::mass, 80.4*GeV, 40*GeV);
      if (imdfound >= 0) {
        const Particle& l = musd_mtfilt[imdfound];
        (l.charge3() > 0 ? _h_Wplus_lepton_eta_mu_dressed : _h_Wminus_lepton_eta_mu_dressed)->fill(l.abseta());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta)
      divide(*_h_Wplus_lepton_eta_el_bare - *_h_Wminus_lepton_eta_el_bare,
             *_h_Wplus_lepton_eta_el_bare + *_h_Wminus_lepton_eta_el_bare,
             _h_W_asym_eta_el_bare);
      divide(*_h_Wplus_lepton_eta_el_dressed - *_h_Wminus_lepton_eta_el_dressed,
             *_h_Wplus_lepton_eta_el_dressed + *_h_Wminus_lepton_eta_el_dressed,
             _h_W_asym_eta_el_dressed);
      divide(*_h_Wplus_lepton_eta_mu_bare - *_h_Wminus_lepton_eta_mu_bare,
             *_h_Wplus_lepton_eta_mu_bare + *_h_Wminus_lepton_eta_mu_bare,
             _h_W_asym_eta_mu_bare);
      divide(*_h_Wplus_lepton_eta_mu_dressed - *_h_Wminus_lepton_eta_mu_dressed,
             *_h_Wplus_lepton_eta_mu_dressed + *_h_Wminus_lepton_eta_mu_dressed,
             _h_W_asym_eta_mu_dressed);

      // Print summary info
      const double xs_pb = crossSection() / picobarn;
      const double sumw = sumOfWeights();
      const double sf = 0.5 * xs_pb / sumw; // 0.5 accounts for rapidity bin width
      scale(_h_Wminus_lepton_eta_el_bare   , sf);
      scale(_h_Wminus_lepton_eta_el_dressed, sf);
      scale(_h_Wminus_lepton_eta_mu_bare   , sf);
      scale(_h_Wminus_lepton_eta_mu_dressed, sf);
      scale(_h_Wplus_lepton_eta_el_bare    , sf);
      scale(_h_Wplus_lepton_eta_el_dressed , sf);
      scale(_h_Wplus_lepton_eta_mu_bare    , sf);
      scale(_h_Wplus_lepton_eta_mu_dressed , sf);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_Wminus_lepton_eta_el_bare;
    Histo1DPtr _h_Wminus_lepton_eta_el_dressed;
    Histo1DPtr _h_Wminus_lepton_eta_mu_bare;
    Histo1DPtr _h_Wminus_lepton_eta_mu_dressed;
    Histo1DPtr _h_Wplus_lepton_eta_el_bare;
    Histo1DPtr _h_Wplus_lepton_eta_el_dressed;
    Histo1DPtr _h_Wplus_lepton_eta_mu_bare;
    Histo1DPtr _h_Wplus_lepton_eta_mu_dressed;
    Estimate1DPtr _h_W_asym_eta_el_bare;
    Estimate1DPtr _h_W_asym_eta_el_dressed;
    Estimate1DPtr _h_W_asym_eta_mu_bare;
    Estimate1DPtr _h_W_asym_eta_mu_dressed;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2011_I928289_W);

}
